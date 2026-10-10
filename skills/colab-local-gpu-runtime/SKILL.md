---
name: colab-local-gpu-runtime
description: Set up, start, verify, and troubleshoot Google Colab local runtimes backed by an NVIDIA GPU on Windows using WSL2, Docker Desktop, and Google's official Colab GPU image. TRIGGER when connecting pdb2reaction, mlmm-toolkit, or another Colab notebook to a Windows-hosted local GPU runtime. SKIP for hosted Colab, Linux-native Jupyter, or ordinary package installation.
---

# Colab Local GPU Runtime

Run the Colab page in the browser while Python and the GPU run locally in Google's official Colab image (Windows with WSL2 and Docker Desktop). `C:\colab-work` is mounted as `/content`, so notebooks and results stay on the disk.

## Steps

1. `powershell -ExecutionPolicy Bypass -File scripts/status.ps1` shows what is missing.
2. First time: `scripts/setup.ps1` installs WSL2 and Docker Desktop, pulls the image, and starts the `colab-gpu` container. Expect a UAC prompt and keep at least 100 GB free.
3. Later: `scripts/start.ps1` starts the container and prints the token URL.
4. In Colab, choose **Connect ▸ Connect to a local runtime**, paste the whole URL including `?token=...`, then run **Installation** and **Launch GUI**.
5. `scripts/stop.ps1` stops the container; `C:\colab-work` is kept.

Run each script with `powershell -ExecutionPolicy Bypass -File`. `setup.ps1` and `start.ps1` copy the notebooks, source ZIPs, and input folders found in the bundle root (two levels above `scripts/`; pass `-BundleRoot` for another layout) into `C:\colab-work`. Keep the notebook's version field at the release tag so Installation installs that release. For an unpublished build, set it to `debug`: Installation then installs from the matching source ZIP next to the notebook, without a file picker.

## Done when

- `status.ps1` shows the same NVIDIA GPU on the host and in the container, `colab-gpu` publishes `127.0.0.1:9000`, and a token URL is printed.
- Every Installation cell finishes, and each Launch GUI cell shows the GUI or emits `application/vnd.jupyter.widget-view+json`.

## When it does not work

- **Connect stays disabled:** paste the full URL with `?token=...`, and open Colab in Chrome or Edge; some embedded browsers block `localhost`, which does not mean the runtime is broken.
- **No GPU in the container:** run `status.ps1` and check the NVIDIA driver on Windows.
- **Docker does not start after the first install:** restart Windows once, open Docker Desktop, and run `start.ps1`.
- **Release not found on PyPI:** use the matching notebook and source ZIP pair and set the version field to `debug`.
- **Files disappear:** save under `/content` (`C:\colab-work`); `drive.mount()` does not work in a local runtime.
- **First check:** the ORB backend needs no login; UMA needs a Hugging Face sign-in and license acceptance.

## Safety

- Run only trusted notebooks; a local runtime can read the files you expose.
- Do not share or commit the token URL.
- Do not delete the image, the container, or `C:\colab-work` unless asked.
