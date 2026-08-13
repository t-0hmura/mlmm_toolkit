---
name: colab-local-gpu-runtime
description: Set up, start, verify, and troubleshoot Google Colab local runtimes backed by an NVIDIA GPU on Windows using WSL2, Docker Desktop, and Google's official Colab GPU image. TRIGGER when connecting pdb2reaction, mlmm-toolkit, or another Colab notebook to a Windows-hosted local GPU runtime. SKIP for hosted Colab, Linux-native Jupyter, or ordinary package installation.
---

# Colab Local GPU Runtime

Run the Colab browser UI against an isolated local NVIDIA GPU container. Map `C:\colab-work` to `/content` so notebooks, source ZIPs, and results persist outside the container.

## Workflow

1. Run `powershell -ExecutionPolicy Bypass -File scripts/status.ps1`.
2. If WSL2, Docker Desktop, the image, or `colab-gpu` is missing, run `powershell -ExecutionPolicy Bypass -File scripts/setup.ps1`. Expect a UAC prompt and reserve at least 100 GB for first-time extraction.
3. On an already configured machine, run `powershell -ExecutionPolicy Bypass -File scripts/start.ps1`.
4. Paste the printed token URL into Colab under **Connect ▸ Connect to a local runtime**.
5. Run **Installation**, then **Launch GUI**.
6. Stop compute with `powershell -ExecutionPolicy Bypass -File scripts/stop.ps1`.

Read `references/colab-connection.md` before connecting or troubleshooting a notebook.

## pdb2reaction and mlmm-toolkit

- Keep the matching notebook and source ZIP together in the distribution root. The scripts copy them to `C:\colab-work` when present.
- For an unpublished bundle, set the notebook version field to `debug`; Installation then consumes the adjacent source ZIP without a file picker.
- Retain the version tag for a published release.
- Use ORB for the easiest no-login smoke test. UMA sign-in and licence acceptance are optional later steps.
- Do not use `drive.mount()` in a local runtime. Exchange files through `/content` / `C:\colab-work`.

## Completion criteria

- Host and container `nvidia-smi` output name the intended GPU.
- `colab-gpu` publishes `127.0.0.1:9000` and `status.ps1` prints a token URL.
- Each Installation cell completes.
- Each Launch GUI cell emits `application/vnd.jupyter.widget-view+json` or visibly renders in Colab.

If an embedded browser blocks `localhost`, test with the user's Chrome or Edge. Do not diagnose the runtime as broken solely from embedded-browser blocking.

## Safety

- Execute only trusted notebooks; local runtime code can access exposed local files.
- Never publish or commit the token URL.
- Never delete the image, container, or `C:\colab-work` unless explicitly requested.
