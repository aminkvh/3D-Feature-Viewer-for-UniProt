# Builds both Chrome (MV3) and Firefox (MV2) extensions in one go.
# Run from the project root:  pwsh ./build-all.ps1

$ErrorActionPreference = 'Stop'
$root = $PSScriptRoot

$files = @(
  'data.js', 'api.js', 'state.js', 'export.js', 'analysis.js', 'algorithms.js',
  'viewer-molstar.js', 'modal.js', 'injector.js', 'content.js',
  'content.css', 'viewer-frame.html', 'viewer-frame.js', 'pocket-worker.js',
  'options.html', 'options.js', 'LICENSE'
)

function Build-Extension($dest, $manifestSrc) {
    if (Test-Path $dest) { Remove-Item $dest -Recurse -Force }
    New-Item -ItemType Directory -Path "$dest/lib"   -Force | Out-Null
    New-Item -ItemType Directory -Path "$dest/icons" -Force | Out-Null
    foreach ($f in $files) { Copy-Item (Join-Path $root $f) (Join-Path $dest $f) -Force }
    Copy-Item (Join-Path $root 'lib/molstar.js')  "$dest/lib/molstar.js"  -Force
    Copy-Item (Join-Path $root 'lib/molstar.css') "$dest/lib/molstar.css" -Force
    Copy-Item (Join-Path $root 'icons/icon48.png')  "$dest/icons/icon48.png"  -Force
    Copy-Item (Join-Path $root 'icons/icon128.png') "$dest/icons/icon128.png" -Force
    Copy-Item (Join-Path $root $manifestSrc) "$dest/manifest.json" -Force
}

Write-Host "Building Chrome (MV3)..." -ForegroundColor Cyan
Build-Extension (Join-Path $root 'chrome-build') 'manifest.json'
Write-Host "  chrome-build/ ready" -ForegroundColor Green

Write-Host "Building Firefox (MV2)..." -ForegroundColor Cyan
Build-Extension (Join-Path $root 'firefox-build') 'manifest.firefox.json'
Write-Host "  firefox-build/ ready" -ForegroundColor Green

Write-Host ""
Write-Host "Chrome:  chrome://extensions  -> Load unpacked -> chrome-build/"
Write-Host "Firefox: about:debugging      -> Load Temporary Add-on -> firefox-build/manifest.json"
