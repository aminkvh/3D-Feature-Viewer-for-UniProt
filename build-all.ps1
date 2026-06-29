# Builds both Chrome (MV3) and Firefox (MV2) extensions and produces release .zip files.
# Run from the project root:  pwsh ./build-all.ps1
#
# Why two manifests:
#   Chrome runs Mol* in a `sandbox` page (a Chrome-only manifest key that allows unsafe-eval).
#   Firefox MV2 doesn't support the sandbox key; Firefox MV3 forbids unsafe-eval on extension
#   pages. So the Firefox build stays on MV2 with a string CSP.

$ErrorActionPreference = 'Stop'
$root = $PSScriptRoot

# Read version from manifest
$version = (Get-Content (Join-Path $root 'manifest.json') | ConvertFrom-Json).version

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

function Zip-Build($srcDir, $zipName) {
    $zipPath = Join-Path $root $zipName
    if (Test-Path $zipPath) { Remove-Item $zipPath -Force }
    Compress-Archive -Path "$srcDir/*" -DestinationPath $zipPath
    Write-Host "  -> $zipName" -ForegroundColor Gray
}

Write-Host "Building Chrome (MV3) v$version..." -ForegroundColor Cyan
$chromeDest = Join-Path $root 'chrome-build'
Build-Extension $chromeDest 'manifest.json'
Zip-Build $chromeDest "chrome-extension-v$version.zip"
Write-Host "  chrome-build/ ready" -ForegroundColor Green

Write-Host "Building Firefox (MV2) v$version..." -ForegroundColor Cyan
$firefoxDest = Join-Path $root 'firefox-build'
Build-Extension $firefoxDest 'manifest.firefox.json'
Zip-Build $firefoxDest "firefox-extension-v$version.zip"
Write-Host "  firefox-build/ ready" -ForegroundColor Green

Write-Host ""
Write-Host "Release packages:" -ForegroundColor Yellow
Write-Host "  chrome-extension-v$version.zip   — Chrome Web Store / load unpacked"
Write-Host "  firefox-extension-v$version.zip  — Mozilla Add-ons / load temporary"
Write-Host ""
Write-Host "Manual load:"
Write-Host "  Chrome:  chrome://extensions  -> Developer mode -> Load unpacked -> chrome-build/"
Write-Host "  Firefox: about:debugging      -> Load Temporary Add-on -> firefox-build/manifest.json"
