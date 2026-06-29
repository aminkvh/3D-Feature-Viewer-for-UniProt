# Builds a Firefox-ready (Manifest V2) copy of the extension into firefox-build/.
#
# Why a separate build: the Chrome version runs Mol* in a `sandbox` page (a Chrome-only
# manifest key) because Mol* needs 'unsafe-eval'. Firefox doesn't support the sandbox key,
# and Firefox MV3 forbids 'unsafe-eval' on extension pages — so the Firefox build is MV2,
# where the page CSP may grant 'unsafe-eval'. The JS source is identical; only the manifest
# differs (manifest.firefox.json -> firefox-build/manifest.json).
#
# Run from the project root:  pwsh ./build-firefox.ps1   (or right-click > Run with PowerShell)

$ErrorActionPreference = 'Stop'
$root = $PSScriptRoot
$dest = Join-Path $root 'firefox-build'

# Clean rebuild
if (Test-Path $dest) { Remove-Item $dest -Recurse -Force }
New-Item -ItemType Directory -Path $dest            | Out-Null
New-Item -ItemType Directory -Path "$dest/lib"      | Out-Null
New-Item -ItemType Directory -Path "$dest/icons"    | Out-Null

# Runtime files the content script, the Mol* iframe, the worker and the options page need.
$files = @(
  'data.js', 'api.js', 'state.js', 'export.js', 'analysis.js', 'algorithms.js',
  'viewer-molstar.js', 'modal.js', 'injector.js', 'content.js',
  'content.css', 'viewer-frame.html', 'viewer-frame.js', 'pocket-worker.js',
  'options.html', 'options.js', 'LICENSE'
)
foreach ($f in $files) { Copy-Item (Join-Path $root $f) (Join-Path $dest $f) -Force }

# Mol* bundle + its CSS
Copy-Item (Join-Path $root 'lib/molstar.js')  "$dest/lib/molstar.js"  -Force
Copy-Item (Join-Path $root 'lib/molstar.css') "$dest/lib/molstar.css" -Force

# Icons referenced by the manifest
Copy-Item (Join-Path $root 'icons/icon48.png')  "$dest/icons/icon48.png"  -Force
Copy-Item (Join-Path $root 'icons/icon128.png') "$dest/icons/icon128.png" -Force

# Firefox MV2 manifest
Copy-Item (Join-Path $root 'manifest.firefox.json') "$dest/manifest.json" -Force

Write-Host "Firefox build ready at: $dest"
Write-Host "Load it via about:debugging > This Firefox > Load Temporary Add-on > pick firefox-build/manifest.json"
