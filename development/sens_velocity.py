import os
import csv
from PIL import Image
import math
import numpy as np

def find_red_pixels(image_path):
    with Image.open(image_path) as img:
        pixels = img.load()
        width, height = img.size
        red_pixels = []
        for x in range(width):
            for y in range(height):
                r, g, b = pixels[x, y]
                if r > 200 and g < 50 and b < 50:
                    red_pixels.append((x, y))
        if red_pixels:
            avg_x = sum(x for x, _ in red_pixels) / len(red_pixels)
            avg_y = sum(y for _, y in red_pixels) / len(red_pixels)
            
            # Round to 0 decimal places (i.e., to integers)
            avg_x = round(avg_x)
            avg_y = round(avg_y)
            
            #print(f"Found red pixels in {image_path}. Average position: ({avg_x}, {avg_y})")
            return (avg_x, avg_y)
    
    print(f"No red pixels found in {image_path}")
    return None

def calculate_velocity(coord1, coord2, pixel_to_mm, frame_rate, frame_diff):
    x1, y1 = coord1
    x2, y2 = coord2
    
    # Full 2D velocity calculation
    pixel_distance = math.sqrt((x2 - x1)**2 + (y2 - y1)**2)
    distance_mm = pixel_distance / pixel_to_mm
    time_elapsed = frame_diff / frame_rate
    velocity_2d_m_s = (distance_mm / 1000) / time_elapsed
    
    # Simplified horizontal velocity calculation
    horizontal_pixel_distance = abs(x2 - x1)
    horizontal_distance_mm = horizontal_pixel_distance / pixel_to_mm
    velocity_horizontal_m_s = (horizontal_distance_mm / 1000) / time_elapsed
    
    # Calculate 5% uncertainty
    uncertainty_2d = round(velocity_2d_m_s * 0.05, 2)
    uncertainty_horizontal = round(velocity_horizontal_m_s * 0.05, 2)
    
    # Round velocities to 2 decimal places
    velocity_2d_m_s = round(velocity_2d_m_s, 2)
    velocity_horizontal_m_s = round(velocity_horizontal_m_s, 2)
    
    print(f"2D Velocity: {velocity_2d_m_s} ± {uncertainty_2d} m/s")
    print(f"Horizontal Velocity: {velocity_horizontal_m_s} ± {uncertainty_horizontal} m/s")
    
    return velocity_2d_m_s, velocity_horizontal_m_s, uncertainty_2d, uncertainty_horizontal

def process_images(folder_path, pixel_to_mm, frame_rate, frame_diff):
    results = {}
    files = os.listdir(folder_path)
    
    # Group files by sensor name
    sensor_files = {}
    for file in files:
        if file.endswith('.jpeg'):
            sensor_name = file.split('_')[0]
            if sensor_name not in sensor_files:
                sensor_files[sensor_name] = []
            sensor_files[sensor_name].append(file)
    
    for sensor_name, sensor_images in sensor_files.items():
        if len(sensor_images) == 2:
            frame1_file = next(f for f in sensor_images if f.endswith('_1.jpeg'))
            frame2_file = next(f for f in sensor_images if f.endswith('_2.jpeg'))
            
            frame1_path = os.path.join(folder_path, frame1_file)
            frame2_path = os.path.join(folder_path, frame2_file)
            
            print(f"Processing sensor: {sensor_name}")
            #print(f"Frame 1 path: {frame1_path}")
            #print(f"Frame 2 path: {frame2_path}")
            
            coord1 = find_red_pixels(frame1_path)
            coord2 = find_red_pixels(frame2_path)
            
            if coord1 and coord2:
                velocity = calculate_velocity(coord1, coord2, pixel_to_mm, frame_rate, frame_diff)
                results[sensor_name] = {
                    'frame1_coord': coord1,
                    'frame2_coord': coord2,
                    'velocity_2d_m_s': velocity[0],
                    'velocity_2d_uncertainty_m_s': velocity[2],
                    'velocity_horizontal_m_s': velocity[1],
                    'velocity_horizontal_uncertainty_m_s': velocity[3],
                    'pixel_to_mm': pixel_to_mm,
                    'frame_rate': frame_rate,
                    'frame_diff': frame_diff
                }
                print(f"Added results for sensor {sensor_name}")
            else:
                print(f"Could not find red pixels in both frames for sensor {sensor_name}")
        else:
            print(f"Incomplete image pair for sensor {sensor_name}")
    
    print(f"Total sensors processed: {len(results)}")
    return results

def calculate_combined_uncertainty(results):
    velocities_2d = [data['velocity_2d_m_s'] for data in results.values()]
    uncertainties_2d = [data['velocity_2d_uncertainty_m_s'] for data in results.values()]
    velocities_horizontal = [data['velocity_horizontal_m_s'] for data in results.values()]
    uncertainties_horizontal = [data['velocity_horizontal_uncertainty_m_s'] for data in results.values()]

    def process_velocity(velocities, uncertainties):
        mean_velocity = np.mean(velocities)
        sd = np.std(velocities, ddof=1)
        propagated_uncertainty = np.sqrt(np.sum(np.array(uncertainties)**2)) / len(uncertainties)
        combined_uncertainty = np.sqrt(sd**2 + propagated_uncertainty**2)
        return mean_velocity, sd, propagated_uncertainty, combined_uncertainty

    results_2d = process_velocity(velocities_2d, uncertainties_2d)
    results_horizontal = process_velocity(velocities_horizontal, uncertainties_horizontal)

    return {
        'mean_velocity_2d': round(results_2d[0], 2),
        'sd_velocity_2d': round(results_2d[1], 2),
        'propagated_uncertainty_2d': round(results_2d[2], 2),
        'combined_uncertainty_2d': round(results_2d[3], 2),
        'mean_velocity_horizontal': round(results_horizontal[0], 2),
        'sd_velocity_horizontal': round(results_horizontal[1], 2),
        'propagated_uncertainty_horizontal': round(results_horizontal[2], 2),
        'combined_uncertainty_horizontal': round(results_horizontal[3], 2)
    }

def save_results_to_csv(results, output_path):
    with open(output_path, 'w', newline='') as csvfile:
        fieldnames = ['sensor_name', 'frame1_x', 'frame1_y', 'frame2_x', 'frame2_y', 
                      'velocity_2d_m_s', 'velocity_2d_uncertainty_m_s', 
                      'velocity_horizontal_m_s', 'velocity_horizontal_uncertainty_m_s', 
                      'pixel_to_mm', 'frame_rate', 'frame_diff']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        
        writer.writeheader()
        for sensor_name, data in results.items():
            writer.writerow({
                'sensor_name': sensor_name,
                'frame1_x': data['frame1_coord'][0],
                'frame1_y': data['frame1_coord'][1],
                'frame2_x': data['frame2_coord'][0],
                'frame2_y': data['frame2_coord'][1],
                'velocity_2d_m_s': data['velocity_2d_m_s'],
                'velocity_2d_uncertainty_m_s': data['velocity_2d_uncertainty_m_s'],
                'velocity_horizontal_m_s': data['velocity_horizontal_m_s'],
                'velocity_horizontal_uncertainty_m_s': data['velocity_horizontal_uncertainty_m_s'],
                'pixel_to_mm': data['pixel_to_mm'],
                'frame_rate': data['frame_rate'],
                'frame_diff': data['frame_diff']
            })
        
        # Add combined uncertainty information
        combined_uncertainty = calculate_combined_uncertainty(results)
        writer.writerow(dict.fromkeys(fieldnames, ''))  # Add an empty row for separation
        for key, value in combined_uncertainty.items():
            writer.writerow({fieldnames[0]: key, fieldnames[1]: value})

        print(f"Wrote {len(results)} rows to CSV, plus combined uncertainty information")

# Main execution
if __name__ == "__main__":
    folder_path = './development/velocity_captures'
    pixel_to_mm = 1.65
    frame_rate = 1000
    frame_diff = 100

    #print(f"Processing images in folder: {folder_path}")
    results = process_images(folder_path, pixel_to_mm, frame_rate, frame_diff)
    output_csv_path = os.path.join(folder_path, 'velocity_results.csv')
    save_results_to_csv(results, output_csv_path)

    print(f"Processing complete. Results saved to {output_csv_path}")

#############
# 
# def process_images(folder_path):
#     print(f"Searching for files in: {folder_path}")
#     if not os.path.exists(folder_path):
#         print(f"Error: Folder {folder_path} does not exist")
#         return
# 
#     all_files = os.listdir(folder_path)
#     print(f"Total files found: {len(all_files)}")
#     
#     for filename in all_files:
#         print(f"Found file: {filename}")
#         if filename.endswith('_1'):
#             sensor_name = filename.split('_')[0]
#             frame1_path = os.path.join(folder_path, filename)
#             frame2_path = os.path.join(folder_path, filename[:-2] + '2')
#             
#             print(f"  Potential sensor: {sensor_name}")
#             print(f"  Frame 1 path: {frame1_path}")
#             print(f"  Frame 2 path: {frame2_path}")
#             
#             if os.path.exists(frame2_path):
#                 print(f"  Matching frame 2 found")
#             else:
#                 print(f"  No matching frame 2 found")
#         else:
#             print(f"  Not a frame 1 file, skipping")
# 
# # Main execution
# if __name__ == "__main__":
#     folder_path = './development/velocity_captures'
#     process_images(folder_path)
