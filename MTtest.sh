#!/bin/bash
# ============================================================
# Automated test script for MyTunesNew (MTconvert / MTgraph)
# ============================================================
# It will:
#   1. Run MTconvert in different modes
#   2. Run MTgraph on the spectrum output
#   3. Verify file existence and sizes
#   4. Produce summary results
# ============================================================
#
# usage: bash MTtest.sh
#
# ============================================================

set -e  # stop on first error
set -o pipefail

WAV="test_440.wav"
SPEC="spectrum.txt"
SPEC_PS="spectrum.ps"
GRAPH_PS="spectrogram.ps"
PIPE_PS="spectrogram_pipe.ps"

echo "=============================================="
echo " MyTunesNew Test Script"
echo "=============================================="
echo

# ---------------------------
# 1. Basic usage check
# ---------------------------
echo "[1] Checking usage message..."
if ./MTconvert >/dev/null 2>&1; then
  echo "❌ Expected usage error, but MTconvert ran silently!"
  exit 1
else
  echo "✅ Usage message printed correctly"
fi
echo

# ---------------------------
# 2. Generate text spectrum only
# ---------------------------
echo "[2] Generating spectrum text only..."
./MTconvert "$WAV" "$SPEC" > run.log 2>&1
if [ -s "$SPEC" ]; then
  echo "✅ Spectrum text created: $(wc -l < "$SPEC") lines"
else
  echo "❌ spectrum.txt missing or empty"
  exit 1
fi
echo

# ---------------------------
# 3. Generate spectrum + PostScript
# ---------------------------
echo "[3] Generating spectrum text + spectrum.ps..."
./MTconvert "$WAV" "$SPEC" "$SPEC_PS" > run_ps.log 2>&1
if [ -s "$SPEC_PS" ]; then
  echo "✅ spectrum.ps created: $(du -h "$SPEC_PS" | cut -f1)"
else
  echo "❌ spectrum.ps missing or empty"
  exit 1
fi
echo

# ---------------------------
# 4. MTgraph standalone test
# ---------------------------
echo "[4] Running MTgraph standalone..."
./MTgraph "$SPEC" "$GRAPH_PS" > graph.log 2>&1
if [ -s "$GRAPH_PS" ]; then
  echo "✅ spectrogram.ps created: $(du -h "$GRAPH_PS" | cut -f1)"
else
  echo "❌ spectrogram.ps missing or empty"
  exit 1
fi
echo

# ---------------------------
# 5. Pipeline test (convert | graph)
# ---------------------------
echo "[5] Testing pipeline mode..."
./MTconvert "$WAV" - | ./MTgraph > "$PIPE_PS" 2> pipe.log
if [ -s "$PIPE_PS" ]; then
  echo "✅ Pipe test successful: $(du -h "$PIPE_PS" | cut -f1)"
else
  echo "❌ Pipe output failed"
  exit 1
fi
echo

# ---------------------------
# 6. Summary
# ---------------------------
echo "=============================================="
echo "✅ All tests passed successfully"
echo "Outputs created:"
ls -lh "$SPEC" "$SPEC_PS" "$GRAPH_PS" "$PIPE_PS"
echo "=============================================="
