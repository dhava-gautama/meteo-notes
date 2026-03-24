# Intel ML Optimization Guide — Smoll Workstation

**Machine**: i7-12700 (12C/20T, P-cores only for compute), 16GB RAM, no GPU
**OS**: Linux 6.17, oneAPI 2025.3
**Date**: 2026-03-24

---

## Python Environment

```bash
# Intel Python — USE THIS for all ML/data science tasks
/opt/intel/oneapi/intelpython/python3.12/bin/python3

# Install packages
/opt/intel/oneapi/intelpython/python3.12/bin/pip install PACKAGE_NAME
```

### Installed Packages
scikit-learn 1.8, sklearnex (oneDAL), XGBoost 3.2, LightGBM 4.6, PyTorch 2.7+cpu, IPEX 2.7+cpu, OpenVINO 2026.0, NNCF, onnx, onnxruntime, xarray, netCDF4, cartopy, matplotlib, pandas, numpy (MKL), scipy, cfgrib

### CPU Features
- **AVX2** — baseline SIMD
- **VNNI** — accelerates INT8 neural network inference
- **No AVX-512, no AMX** — BFloat16 weight prepacking not supported

---

## Required Code Patterns

### scikit-learn (6x speedup)
```python
from sklearnex import patch_sklearn
patch_sklearn()  # MUST be before any sklearn import

from sklearn.ensemble import RandomForestClassifier
# ... use sklearn normally, oneDAL accelerates transparently
```

### PyTorch Training (10-50% speedup)
```python
import torch
import intel_extension_for_pytorch as ipex

model = MyModel()
optimizer = torch.optim.Adam(model.parameters())
model, optimizer = ipex.optimize(model, optimizer=optimizer)
# ... train normally
```

### PyTorch Inference with OpenVINO INT8 (7-10x speedup)
```python
import torch
import openvino as ov
import nncf

# 1. Export trained model to ONNX
model.eval()
dummy = torch.randn(1, ...)  # match your input shape
torch.onnx.export(model, dummy, "model.onnx", opset_version=17)

# 2. Load into OpenVINO
core = ov.Core()
ov_model = core.read_model("model.onnx")

# 3. Quantize to INT8 (needs calibration data)
def transform_fn(data_item):
    return {0: data_item}

calibration_data = [np.random.randn(...).astype(np.float32) for _ in range(100)]
quantized = nncf.quantize(
    ov_model,
    nncf.Dataset(calibration_data, transform_fn),
    preset=nncf.QuantizationPreset.PERFORMANCE,
    fast_bias_correction=True
)

# 4. Compile and run
compiled = core.compile_model(quantized, "CPU")
infer_req = compiled.create_infer_request()
result = infer_req.infer({0: input_numpy_array})

# Optional: save compiled model for reuse
ov.save_model(quantized, "model_int8.xml")
```

---

## Benchmark Results

### scikit-learn: Vanilla vs Intel (sklearnex/oneDAL)
Dataset: 50,000 samples x 30 features

| Task | Vanilla | Intel | Speedup |
|------|---------|-------|---------|
| RandomForest (200 trees) | 2.98s | 0.36s | **8.3x** |
| RandomForest Regressor | 23.24s | 3.75s | **6.2x** |
| Logistic Regression | 0.080s | 0.024s | **3.4x** |
| KNN (k=5) | 0.445s | 0.073s | **6.1x** |
| KMeans (k=8) | 0.394s | 0.219s | 1.8x |
| PCA (10 components) | 0.014s | 0.007s | 2.0x |
| **TOTAL** | **27.17s** | **4.44s** | **6.1x** |

### Generic DL Models: PyTorch vs IPEX vs OpenVINO

| Model | PyTorch | IPEX | OV FP32 | OV INT8 | Speedup |
|-------|---------|------|---------|---------|---------|
| U-Net 64x64 | 0.271s | 0.292s | 0.104s | **0.030s** | **9.0x** |
| U-Net 128x128 | 0.530s | 0.442s | 0.217s | **0.078s** | **6.8x** |
| U-Net 256x256 | 1.174s | 0.758s | 0.436s | **0.177s** | **6.6x** |
| MLP 4k x 50 | 0.015s | 0.014s | 0.007s | **0.002s** | **7.2x** |
| MLP 65k x 50 | 0.388s | 0.385s | 0.129s | **0.040s** | **9.6x** |

### Weather Architectures: PyTorch vs OpenVINO INT8
Input: batch=8, 10 timesteps, 1 channel, 64x64

| Model | Params | PyTorch | OV INT8 | Speedup |
|-------|--------|---------|---------|---------|
| ConvLSTM-ED (64 hidden) | 300K | 1.051s | 0.414s | **2.5x** |
| BiEF (64 hidden) | 467K | 1.370s | 0.628s | **2.2x** |
| VBiEF (64h/32z, beta=0.001) | 465K | 1.377s | 0.641s | **2.1x** |
| CNN-GRU-MS (256 hidden) | 1.71M | 0.026s | 0.004s | **6.2x** |

At 128x128:

| Model | PyTorch | OV INT8 | Speedup |
|-------|---------|---------|---------|
| ConvLSTM-ED 128x128 | 1.919s | 0.850s | **2.3x** |
| BiEF 128x128 | 2.879s | 1.282s | **2.2x** |
| CNN-GRU-MS 128x128 | 0.029s | 0.005s | **6.5x** |

### Advanced Architectures: SmaAt-UNet, ViT, Swin, AFNO
Input: batch=8, 10 channels, 64x64

| Model | Params | PyTorch | IPEX | OV FP32 | OV INT8 | Speedup |
|-------|--------|---------|------|---------|---------|---------|
| SmaAt-UNet (CBAM) | 7.72M | 0.141s | 0.227s | 0.054s | **0.017s** | **8.5x** |
| Swin Transformer (192d, 12L) | 5.37M | 0.294s | 0.318s | 0.048s | **0.027s** | **10.9x** |
| ViT (256d, 6L) | 4.94M | 0.023s | 0.025s | 0.010s | **0.004s** | **5.2x** |
| AFNO/FourCastNet (256d, 6L) | 3.56M | 0.024s | 0.024s | FAIL | FAIL | 1.0x |

At 128x128 (batch=4):

| Model | PyTorch | IPEX | OV FP32 | OV INT8 | Speedup |
|-------|---------|------|---------|---------|---------|
| SmaAt-UNet 128x128 | 0.284s | 0.318s | 0.116s | **0.036s** | **7.8x** |
| Swin Transformer 128x128 | 1.926s | 1.906s | 0.192s | **0.104s** | **18.6x** |
| ViT 128x128 | 0.064s | 0.065s | 0.022s | **0.010s** | **6.2x** |
| AFNO 128x128 | 0.043s | 0.042s | FAIL | FAIL | 1.0x |

Notes:
- **Swin Transformer gets the biggest OpenVINO boost** (18.6x at 128x128) — windowed attention is perfectly suited for INT8/VNNI
- **AFNO cannot use OpenVINO** — `fft_rfft2` is not supported in ONNX export. Stuck at PyTorch speed.
- **IPEX actually hurts** SmaAt-UNet and Swin (slower than vanilla PyTorch)

### WRF 4.7.1 Core Scaling
Domain: 100x100, dx=9km, 45 levels, 6h forecast

| Cores | Wall Time | Speedup | Efficiency |
|-------|-----------|---------|------------|
| 1 | 657s | 1.00x | 100% |
| 2 | 396s | 1.66x | 83% |
| 4 | 294s | 2.23x | 56% |
| 6 | 255s | 2.58x | 43% |
| 8 | 236s | 2.78x | 35% |

Sweet spot: **6 cores** for 100x100 domains. Use 8 for larger domains (300x300+).

---

## Decision Guide

| Task | Tool | Expected Speedup |
|------|------|-----------------|
| Classical ML (RF, KNN, SVM, KMeans) | sklearnex | 6-8x |
| Gradient boosting | XGBoost/LightGBM (native) | MKL-accelerated |
| DL training | PyTorch + IPEX | 10-50% |
| DL inference (CNN, MLP) | OpenVINO INT8 | **7-10x** |
| DL inference (recurrent/ConvLSTM) | OpenVINO INT8 | 2-2.5x |
| DL inference (CNN-GRU hybrid) | OpenVINO INT8 | **6x** |
| DL inference (Swin Transformer) | OpenVINO INT8 | **11-19x** |
| DL inference (SmaAt-UNet/attention) | OpenVINO INT8 | **8-9x** |
| DL inference (ViT) | OpenVINO INT8 | **5-6x** |
| DL inference (AFNO/spectral) | PyTorch only | no OV support (FFT ops) |
| NumPy/SciPy heavy compute | Already MKL-accelerated | automatic |

### Architecture Recommendation for CPU Deployment

**Tier 1 — Best on CPU (OpenVINO INT8)**
- **Swin Transformer**: 18.6x speedup at 128x128 — windowed attention is ideal for VNNI
- **SmaAt-UNet**: 8.5x — proven for precipitation nowcasting, CNN+CBAM optimizes well
- **CNN-GRU-MS**: 6.2x and 50x faster than ConvLSTM — best recurrent architecture

**Tier 2 — Good on CPU**
- **ViT**: 5-6x speedup, simpler than Swin but less optimizable
- **U-Net (plain)**: 7-10x — baseline workhorse, very reliable

**Tier 3 — Poor on CPU**
- **ConvLSTM variants** (ED, BiEF, VBiEF): Only 2-2.5x — sequential loops limit optimization
- **AFNO/FourCastNet**: No OpenVINO support (FFT not in ONNX) — can't be accelerated beyond PyTorch

**Key insight**: For this CPU (i7-12700 + VNNI), prefer architectures with dense matrix ops (convolutions, linear layers, attention) over sequential/spectral ops. Swin > SmaAt-UNet > ViT > CNN-GRU > ConvLSTM > AFNO for deployment speed.

### Workflow
```
Train (PyTorch + IPEX) → Export (ONNX) → Quantize (NNCF INT8) → Deploy (OpenVINO)
```

---

## Limitations
- **16GB RAM**: WRF compilation must use `-j 1` (RSL_LITE files need 7+ GB each at -O3)
- **No GPU**: All inference is CPU-only
- **No BFloat16**: i7-12700 lacks AVX-512/AMX, IPEX weights_prepack fails
- **ConvLSTM ONNX export**: Works but sequential loops can't be fully optimized by OpenVINO