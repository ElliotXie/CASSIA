#!/usr/bin/env python3
"""
Image Optimization Script for CASSIA Documentation

USAGE:
  1. Drop your PNG/JPG images into public/images/
  2. Run: python scripts/optimize-images.py
  3. Done! Images are optimized and originals backed up.

The script will:
  - Detect new unoptimized images (PNG/JPG)
  - Back up originals to public/images_original/
  - Convert to WebP (resized if needed)
  - Keep whichever is smaller (WebP or original)
  - Clean up automatically
"""

import os
import shutil
from pathlib import Path
from PIL import Image

# Configuration
MAX_WIDTH = 1400  # Max width in pixels (supports retina displays)
WEBP_QUALITY = 92  # High quality, visually lossless

# Directories
SCRIPT_DIR = Path(__file__).parent
IMAGES_DIR = SCRIPT_DIR.parent / "public" / "images"
BACKUP_DIR = SCRIPT_DIR.parent / "public" / "images_original"

# Supported formats
IMAGE_EXTENSIONS = {'.png', '.jpg', '.jpeg'}


def get_file_size_kb(path):
    """Get file size in KB."""
    return os.path.getsize(path) / 1024


def optimize_image(input_path, output_path, max_width=MAX_WIDTH, quality=WEBP_QUALITY):
    """
    Optimize a single image:
    - Resize if wider than max_width (maintaining aspect ratio)
    - Convert to WebP format
    """
    with Image.open(input_path) as img:
        original_size = img.size

        # Handle different color modes
        if img.mode in ('RGBA', 'LA', 'P'):
            if img.mode == 'P':
                img = img.convert('RGBA')
        elif img.mode != 'RGB':
            img = img.convert('RGB')

        # Resize if needed
        if img.width > max_width:
            ratio = max_width / img.width
            new_height = int(img.height * ratio)
            img = img.resize((max_width, new_height), Image.Resampling.LANCZOS)

        # Save as WebP
        img.save(output_path, 'WEBP', quality=quality, method=6)

        return original_size, img.size


def find_new_images():
    """Find PNG/JPG images that haven't been optimized yet."""
    new_images = []

    for img_path in IMAGES_DIR.iterdir():
        if img_path.suffix.lower() not in IMAGE_EXTENSIONS:
            continue

        # Check if this image is already backed up (meaning it's been processed)
        backup_path = BACKUP_DIR / img_path.name
        if backup_path.exists():
            continue

        new_images.append(img_path)

    return new_images


def process_image(img_path):
    """
    Process a single image:
    1. Back up original
    2. Create WebP version
    3. Compare sizes
    4. Keep smaller version, delete other

    Returns: dict with results
    """
    original_kb = get_file_size_kb(img_path)

    # Ensure backup directory exists
    BACKUP_DIR.mkdir(parents=True, exist_ok=True)

    # Back up original
    backup_path = BACKUP_DIR / img_path.name
    shutil.copy2(img_path, backup_path)

    # Create WebP version
    webp_name = img_path.stem + '.webp'
    webp_path = IMAGES_DIR / webp_name

    try:
        orig_dims, new_dims = optimize_image(img_path, webp_path)
        webp_kb = get_file_size_kb(webp_path)

        # Decide which to keep
        if webp_kb < original_kb:
            # WebP is smaller - delete original, keep WebP
            img_path.unlink()
            kept_format = 'webp'
            kept_kb = webp_kb
            final_name = webp_name
        else:
            # Original is smaller - delete WebP, keep original
            webp_path.unlink()
            kept_format = img_path.suffix[1:]  # Remove the dot
            kept_kb = original_kb
            final_name = img_path.name

        savings = ((original_kb - kept_kb) / original_kb) * 100

        return {
            'success': True,
            'original_name': img_path.name,
            'final_name': final_name,
            'original_kb': original_kb,
            'kept_kb': kept_kb,
            'kept_format': kept_format,
            'savings': savings,
            'dimensions': f"{orig_dims[0]}x{orig_dims[1]} -> {new_dims[0]}x{new_dims[1]}"
        }

    except Exception as e:
        return {
            'success': False,
            'original_name': img_path.name,
            'error': str(e)
        }


def main():
    print("=" * 65)
    print("  CASSIA Documentation - Image Optimizer")
    print("=" * 65)
    print()
    print(f"  Settings: max {MAX_WIDTH}px width, WebP quality {WEBP_QUALITY}")
    print(f"  Images folder: {IMAGES_DIR}")
    print(f"  Backup folder: {BACKUP_DIR}")
    print()

    # Find new images to process
    new_images = find_new_images()

    if not new_images:
        print("  No new images to optimize!")
        print()
        print("  To add new images:")
        print("    1. Drop PNG/JPG files into public/images/")
        print("    2. Run this script again")
        print()

        # Show current status
        current_images = list(IMAGES_DIR.iterdir())
        if current_images:
            total_size = sum(f.stat().st_size for f in current_images if f.is_file())
            print(f"  Current images: {len([f for f in current_images if f.is_file()])} files, {total_size/1024:.0f} KB total")
        return

    print(f"  Found {len(new_images)} new image(s) to optimize:")
    print()

    # Process each image
    results = []
    for img_path in sorted(new_images):
        print(f"  Processing: {img_path.name}...", end=" ", flush=True)
        result = process_image(img_path)
        results.append(result)

        if result['success']:
            print(f"OK -> {result['final_name']}")
        else:
            print(f"ERROR: {result['error']}")

    # Summary
    print()
    print("-" * 65)
    print(f"  {'Image':<30} {'Original':>10} {'Final':>10} {'Saved':>10}")
    print("-" * 65)

    total_original = 0
    total_final = 0

    for r in results:
        if r['success']:
            total_original += r['original_kb']
            total_final += r['kept_kb']
            print(f"  {r['final_name']:<30} {r['original_kb']:>8.0f}KB {r['kept_kb']:>8.0f}KB {r['savings']:>8.1f}%")
        else:
            print(f"  {r['original_name']:<30} {'ERROR':>10}")

    if total_original > 0:
        total_savings = ((total_original - total_final) / total_original) * 100
        print("-" * 65)
        print(f"  {'TOTAL':<30} {total_original:>8.0f}KB {total_final:>8.0f}KB {total_savings:>8.1f}%")

    print()
    print("  Done! Originals backed up to: images_original/")
    print()
    print("  IMPORTANT: Update your markdown to use the new filename!")
    print("  Example: ![Description](/images/your-image.webp)")
    print()


if __name__ == "__main__":
    main()
