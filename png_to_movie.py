# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Thu Nov 28 23:16:44 2024
"""

from PIL import Image, ImageDraw, ImageFont
import cv2
import os

def pngs_to_animation(image_folder, output_file, selected_files=None, duration=500, frame_text=False, output_type="gif"):
    """
    Converts selected PNG images in a folder into a single GIF or MP4 video with optional text on each frame.

    Args:
    - image_folder (str): Path to the folder containing PNG images.
    - output_file (str): Path to save the output file (GIF or MP4).
    - selected_files (list): List of specific image filenames to include. If None, all PNGs in the folder are used.
    - duration (int): Frame duration in milliseconds for GIF or FPS for MP4.
    - frame_text (bool): If True, add frame numbers as text to each frame.
    - output_type (str): Either "gif" for GIF output or "video" for MP4 video output.
    """
    # Check if folder exists
    if not os.path.exists(image_folder):
        raise FileNotFoundError(f"The folder '{image_folder}' does not exist.")

    all_images = [img for img in os.listdir(image_folder) if img.endswith(".png")]
    all_images.sort()

    if selected_files is not None:
        images = [img for img in selected_files if img in all_images]
    else:
        images = all_images

    if not images:
        raise ValueError("No images found or selected to create the animation.")

    if output_type.lower() == "gif":
        # Create GIF
        frames = []
        for i, img_name in enumerate(images):
            img_path = os.path.join(image_folder, img_name)
            img = Image.open(img_path).convert("RGBA")  # Ensure RGBA mode

            # if frame_text:
            #     draw = ImageDraw.Draw(img)
            #     font = ImageFont.load_default()  # Use default font
            #     draw.text((10, 10), f"Frame: {i}", font=font, fill="black")

            frames.append(img)

        frames[0].save(
            output_file,
            save_all=True,
            append_images=frames[1:],
            duration=duration,
            loop=0,
            optimize=True
        )
        print(f"GIF saved at {output_file}")

    elif output_type.lower() == "video":
        # Create MP4 video
        first_frame = cv2.imread(os.path.join(image_folder, images[0]))
        height, width, _ = first_frame.shape
        fps = 1000 // duration  # Convert duration to frames per second

        video_writer = cv2.VideoWriter(
            output_file, cv2.VideoWriter_fourcc(*"mp4v"), fps, (width, height)
        )

        for i, img_name in enumerate(images):
            img_path = os.path.join(image_folder, img_name)
            frame = cv2.imread(img_path)

            if frame_text:
                cv2.putText(
                    frame, f"Frame: {i}",
                    (10, 30), cv2.FONT_HERSHEY_SIMPLEX,
                    1, (0, 0, 255), 2, cv2.LINE_AA
                )

            video_writer.write(frame)

        video_writer.release()
        print(f"Video saved at {output_file}")

    else:
        raise ValueError("Invalid output_type. Choose 'gif' or 'video'.")




# Example usage
inpath=r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\Classes\Math_Fun\winding_transformation\Images'

selected_images = [f'winding_{i}.png' for i in range(0,1000,1)]
pngs_to_animation(inpath, os.path.join(inpath, "outputt.mp4"),
                  selected_files=selected_images,
                  duration=100,
                  output_type="video")

