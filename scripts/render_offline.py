#!/usr/bin/env python3

# This script is used to render the simulation output into a video.
import argparse
import os
import glob
import subprocess
from typing import Tuple

import numpy as np


# This is used to import the Image module from the PIL library. This one has been an issue for some reason which is why it has its own try/except block.
try:
    from PIL import Image
except Exception as e:
    Image = None


# This function is used to load the simulation output into a numpy array.
def load_slice(npy_path: str) -> np.ndarray:
    arr = np.load(npy_path)
    # Normalize to [0, 255] for visualization
    vmin, vmax = np.min(arr), np.max(arr)
    if vmax <= vmin:
        vmax = vmin + 1e-6
    norm = (arr - vmin) / (vmax - vmin)
    img = (norm * 255.0).astype(np.uint8)
    return img


# This function is used to write the frame to a png file.
def write_frame_png(img: np.ndarray, out_path: str, scale: int = 2) -> None:
    if Image is None:
        raise RuntimeError("Pillow is required to write PNGs (pip install Pillow)")
    pil = Image.fromarray(img, mode="L")
    if scale != 1:
        pil = pil.resize((pil.width * scale, pil.height * scale), resample=Image.NEAREST)
    pil.save(out_path)


# This function is used to encode the frames into a video using ffmpeg.
def ffmpeg_encode(frames_glob: str, fps: int, out_video: str) -> None:
    # Requires ffmpeg in PATH
    cmd = [
        "ffmpeg", "-y",
        "-framerate", str(fps),
        "-pattern_type", "glob",
        "-i", frames_glob,
        "-c:v", "libx264", "-pix_fmt", "yuv420p",
        out_video,
    ]
    subprocess.run(cmd, check=True)


# This function is used to collect the steps from the simulation output.
def collect_steps(outdir: str) -> list:
    steps = sorted(glob.glob(os.path.join(outdir, "step_*")))
    return steps


# This function is used to render the sequence of frames into a video.
def render_sequence(input_dir: str, theta: int, field: str, fps: int, scale: int, output: str) -> None:
    steps = collect_steps(input_dir)
    if not steps:
        raise RuntimeError(f"No steps found under {input_dir}. Expected subdirs like step_000000/")

    frames_dir = os.path.join(input_dir, f".frames_th{theta}_{field}")
    os.makedirs(frames_dir, exist_ok=True)

    for idx, step_dir in enumerate(steps):
        npy = os.path.join(step_dir, f"th{theta}_{field}.npy")
        if not os.path.exists(npy):
            # Try legacy format for backward compatibility
            legacy_npy = os.path.join(step_dir, f"tracer_slice_th{theta}.npy")
            if os.path.exists(legacy_npy):
                npy = legacy_npy
            else:
                print(f"Warning: Missing {npy}, skipping frame {idx}")
                continue
        img = load_slice(npy)
        frame_path = os.path.join(frames_dir, f"frame_{idx:06d}.png")
        write_frame_png(img, frame_path, scale=scale)

    frames_glob = os.path.join(frames_dir, "frame_*.png")
    ffmpeg_encode(frames_glob, fps, output)


# This function is used to parse the command line arguments and render the sequence of frames into a video.
def main() -> None:
    p = argparse.ArgumentParser(description="Offline renderer for exported theta-slice sequences")
    p.add_argument("--input", required=True, help="Export directory (e.g., data/exports)")
    p.add_argument("--theta", type=int, default=0, help="Theta index to render (default: 0)")
    p.add_argument("--field", default="tracer", help="Field to render (tracer, theta, qv, qc, qr, qh, qg, u, v, w, rho, p, radar)")
    p.add_argument("--fps", type=int, default=30, help="Output video FPS")
    p.add_argument("--scale", type=int, default=2, help="Integer pixel scale for frames")
    p.add_argument("--output", required=True, help="Output video file (e.g., out.mp4)")
    args = p.parse_args()

    render_sequence(args.input, args.theta, args.field, args.fps, args.scale, args.output)


if __name__ == "__main__":
    main()


