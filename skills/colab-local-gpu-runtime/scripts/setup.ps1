[CmdletBinding()]
param(
    [string]$WorkDir = 'C:\colab-work',
    [string]$ContainerName = 'colab-gpu',
    [int]$HostPort = 9000,
    [string]$BundleRoot = ''
)

$ErrorActionPreference = 'Stop'
if (-not $BundleRoot) { $BundleRoot = (Resolve-Path (Join-Path $PSScriptRoot '..\..')).Path }

function Test-Administrator {
    $identity = [Security.Principal.WindowsIdentity]::GetCurrent()
    $principal = [Security.Principal.WindowsPrincipal]::new($identity)
    $principal.IsInRole([Security.Principal.WindowsBuiltInRole]::Administrator)
}

if (-not (Test-Administrator)) {
    $arguments = @('-NoProfile', '-ExecutionPolicy', 'Bypass', '-File', ('"' + $PSCommandPath + '"'), '-WorkDir', ('"' + $WorkDir + '"'), '-ContainerName', $ContainerName, '-HostPort', $HostPort, '-BundleRoot', ('"' + $BundleRoot + '"'))
    $process = Start-Process powershell.exe -Verb RunAs -ArgumentList $arguments -WindowStyle Hidden -Wait -PassThru
    exit $process.ExitCode
}

if (-not (Get-Command nvidia-smi.exe -ErrorAction SilentlyContinue)) {
    throw 'NVIDIA driver not found. Install the current NVIDIA driver, reboot, and rerun.'
}

Write-Host '[1/6] Enabling WSL2 Windows features...'
foreach ($feature in @('Microsoft-Windows-Subsystem-Linux', 'VirtualMachinePlatform')) {
    & dism.exe /online /enable-feature "/featurename:$feature" /all /norestart | Out-Host
    if ($LASTEXITCODE -notin @(0, 3010)) { throw "Failed to enable $feature (exit $LASTEXITCODE)." }
}

Write-Host '[2/6] Installing WSL and Docker Desktop when missing...'
$wingetCommand = Get-Command winget.exe -ErrorAction SilentlyContinue
$winget = if ($wingetCommand) { $wingetCommand.Source } else { Join-Path $env:LOCALAPPDATA 'Microsoft\WindowsApps\winget.exe' }
if (-not (Test-Path -LiteralPath $winget)) { throw 'winget is unavailable. Install Microsoft App Installer and rerun.' }
if (-not (& wsl.exe --version 2>$null)) {
    & $winget install --exact --id Microsoft.WSL --accept-source-agreements --accept-package-agreements --silent
    if ($LASTEXITCODE -ne 0) { throw "WSL installation failed (exit $LASTEXITCODE)." }
}
$dockerExe = 'C:\Program Files\Docker\Docker\resources\bin\docker.exe'
if (-not (Test-Path -LiteralPath $dockerExe)) {
    & $winget install --exact --id Docker.DockerDesktop --accept-source-agreements --accept-package-agreements --silent
    if ($LASTEXITCODE -ne 0) { throw "Docker Desktop installation failed (exit $LASTEXITCODE)." }
}

Write-Host '[3/6] Starting Docker Desktop...'
$dockerDesktop = 'C:\Program Files\Docker\Docker\Docker Desktop.exe'
if (-not (Get-Process 'Docker Desktop' -ErrorAction SilentlyContinue)) { Start-Process -FilePath $dockerDesktop -WindowStyle Hidden }
$deadline = (Get-Date).AddMinutes(3)
do {
    $previous = $ErrorActionPreference; $ErrorActionPreference = 'SilentlyContinue'
    & $dockerExe info --format '{{.ServerVersion}}' *> $null; $dockerReady = $LASTEXITCODE
    $ErrorActionPreference = $previous
    if ($dockerReady -eq 0) { break }
    Start-Sleep -Seconds 5
} while ((Get-Date) -lt $deadline)
if ($dockerReady -ne 0) { throw 'Docker did not start. Restart Windows once, open Docker Desktop, and rerun.' }

Write-Host '[4/6] Pulling the official Colab GPU image...'
$image = 'asia-docker.pkg.dev/colab-images/public/runtime'
if (-not (& $dockerExe image ls $image --format '{{.ID}}')) {
    $driveName = [IO.Path]::GetPathRoot($WorkDir).TrimEnd('\').TrimEnd(':')
    if ((Get-PSDrive -Name $driveName).Free -lt 100GB) { throw 'First extraction needs at least 100 GB free.' }
}
& $dockerExe pull $image
if ($LASTEXITCODE -ne 0) { throw "Image pull failed (exit $LASTEXITCODE)." }

Write-Host '[5/6] Preparing persistent files...'
New-Item -ItemType Directory -Force -Path $WorkDir | Out-Null
foreach ($name in @('pdb2reaction_colab.ipynb','mlmm_colab.ipynb','pdb2reaction-src.zip','mlmm-toolkit-src.zip','pdb2reaction_inputs','mlmm_inputs')) {
    $source = Join-Path $BundleRoot $name
    if (Test-Path -LiteralPath $source) { Copy-Item -LiteralPath $source -Destination $WorkDir -Recurse -Force }
}

Write-Host '[6/6] Starting the runtime...'
$existing = & $dockerExe ps -a --format '{{.Names}}' | Where-Object { $_ -eq $ContainerName }
if ($existing) { & $dockerExe start $ContainerName | Out-Null }
else { & $dockerExe run -d --name $ContainerName --gpus all --restart unless-stopped -p "127.0.0.1:${HostPort}:8080" -v "${WorkDir}:/content" $image | Out-Null }
Start-Sleep -Seconds 8
& powershell.exe -NoProfile -ExecutionPolicy Bypass -File (Join-Path $PSScriptRoot 'status.ps1') -ContainerName $ContainerName -HostPort $HostPort
