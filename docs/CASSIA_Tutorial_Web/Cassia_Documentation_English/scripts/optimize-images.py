#!/usr/bin/env python3
"""
Image Optimization Script for CASSIA Documentation
Converts images to WebP format with resizing for optimal web performance.
"""

import os
from pathlib import Path
from PIL import Image

# Configuration
MAX_WIDTH = 1400  # Max width in pixels (supports retina displays)
WEBP_QUALITY = 92  # High quality, visually lossless
IMAGES_DIR = Path(__file__).parent.parent / "public" / "images"
BACKUP_DIR = Path(__file__).parent.parent / "public" / "images_original"

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

        # Convert RGBA to RGB if necessary (WebP supports both, but RGB is smaller)
        if img.mode in ('RGBA', 'LA', 'P'):
            # Keep alpha channel for transparency
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

def main():
    print("=" * 60)
    print("CASSIA Documentation Image Optimization")
    print("=" * 60)
    print(f"\nSettings:")
    print(f"  Max width: {MAX_WIDTH}px")
    print(f"  WebP quality: {WEBP_QUALITY}")
    print(f"  Source: {BACKUP_DIR}")
    print(f"  Output: {IMAGES_DIR}")
    print()

    # Get all images from backup directory
    image_extensions = {'.png', '.jpg', '.jpeg'}
    images = [f for f in BACKUP_DIR.iterdir() if f.suffix.lower() in image_extensions]

    if not images:
        print("No images found in backup directory!")
        return

    total_original = 0
    total_optimized = 0
    results = []

    print(f"Processing {len(images)} images...\n")
    print(f"{'Image':<40} {'Original':>10} {'Optimized':>10} {'Savings':>10}")
    print("-" * 70)

    for img_path in sorted(images):
        # Output path with .webp extension
        output_name = img_path.stem + '.webp'
        output_path = IMAGES_DIR / output_name

        # Get original size
        original_kb = get_file_size_kb(img_path)
        total_original += original_kb

        # Optimize
        try:
            orig_dims, new_dims = optimize_image(img_path, output_path)
            optimized_kb = get_file_size_kb(output_path)
            total_optimized += optimized_kb

            savings = ((original_kb - optimized_kb) / original_kb) * 100

            print(f"{img_path.name:<40} {original_kb:>8.0f}KB {optimized_kb:>8.0f}KB {savings:>8.1f}%")

            results.append({
                'name': img_path.name,
                'original_kb': original_kb,
                'optimized_kb': optimized_kb,
                'savings': savings,
                'orig_dims': orig_dims,
                'new_dims': new_dims
            })

        except Exception as e:
            print(f"{img_path.name:<40} ERROR: {e}")

    # Summary
    print("-" * 70)
    total_savings = ((total_original - total_optimized) / total_original) * 100
    print(f"{'TOTAL':<40} {total_original:>8.0f}KB {total_optimized:>8.0f}KB {total_savings:>8.1f}%")
    print()
    print(f"Total saved: {total_original - total_optimized:.0f} KB ({total_savings:.1f}%)")
    print()
    print("Optimization complete!")
    print(f"Original images preserved in: {BACKUP_DIR}")

if __name__ == "__main__":
    main()
