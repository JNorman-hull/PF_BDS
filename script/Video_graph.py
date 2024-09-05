import cv2
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_agg import FigureCanvasAgg
from scipy import interpolate

def load_video(video_path):
    return cv2.VideoCapture(video_path)

def load_csv_data(csv_path):
    return pd.read_csv(csv_path)

def interpolate_data(data, original_rate, target_rate):
    original_time = np.arange(len(data)) / original_rate
    target_time = np.arange(0, len(data) / original_rate, 1 / target_rate)
    f = interpolate.interp1d(original_time, data, kind='cubic', fill_value='extrapolate')
    return f(target_time)

def create_graphs(pressure_data, accmag_data, current_index, data_hz, window_size=0.3):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 2), dpi=300)  # Reverted DPI to 100
    
    time = np.arange(len(pressure_data)) / data_hz
    current_time = current_index / data_hz
    
    # Pressure graph
    ax1.plot(time, pressure_data, 'r-', linewidth=1)  # Reverted line width to 1
    ax1.set_xlim(current_time - window_size/2, current_time + window_size/2)
    pres_min, pres_max = pressure_data.min(), pressure_data.max()
    y_range = pres_max - pres_min
    ax1.set_ylim(pres_min - 0.05*y_range, pres_max + 0.05*y_range)
    ax1.set_ylabel('Pressure (mbar)', fontsize=10)  # Increased font size
    ax1.set_xticks([])
    ax1.tick_params(axis='y', labelsize=8)  # Increased tick label size
    ax1.axvline(x=current_time, color='g', linestyle='--', linewidth=1)
    current_pres = pressure_data[current_index]
    ax1.text(current_time, current_pres, f'{current_pres:.2f}', 
             verticalalignment='center', horizontalalignment='left', fontsize=12)
    
    # Acceleration magnitude graph
    ax2.plot(time, accmag_data, 'b-', linewidth=1)  # Reverted line width to 1
    ax2.set_xlim(current_time - window_size/2, current_time + window_size/2)
    acc_min, acc_max = accmag_data.min(), accmag_data.max()
    y_range = acc_max - acc_min
    ax2.set_ylim(acc_min - 0.05*y_range, acc_max + 0.05*y_range)
    ax2.set_ylabel('Acceleration magnitude\n(m/s\u00b2)', fontsize=10)  # Increased font size
    ax2.set_xticks([])
    ax2.tick_params(axis='y', labelsize=8)  # Increased tick label size
    ax2.axvline(x=current_time, color='g', linestyle='--', linewidth=1)
    current_acc = accmag_data[current_index]
    ax2.text(current_time, current_acc, f'{current_acc:.2f}', 
             verticalalignment='center', horizontalalignment='left', fontsize=12)
    
    plt.tight_layout()
    
    canvas = FigureCanvasAgg(fig)
    canvas.draw()
    graph = np.frombuffer(canvas.buffer_rgba(), dtype='uint8')
    graph = graph.reshape(canvas.get_width_height()[::-1] + (4,))[:,:,:3]
    plt.close(fig)
    return cv2.cvtColor(graph, cv2.COLOR_RGB2BGR)
  

def overlay_graph(frame, graph):
    h, w, _ = frame.shape
    graph_resized = cv2.resize(graph, (w, h//4))
    frame[-h//4:, :] = graph_resized
    return frame

def adjust_frame(frame, video_height_adjust):
    h, w, _ = frame.shape
    adjusted_frame = np.zeros((h, w, 3), dtype=np.uint8)
    adjusted_frame[:-video_height_adjust] = frame[video_height_adjust:]
    return adjusted_frame

def process_video(video_path, csv_path, output_path, nadir_frame, nadir_index, data_offset=0, video_height_adjust=0, video_text=None):
    if video_text is None:
        video_text = []
    
    video = load_video(video_path)
    data = load_csv_data(csv_path)
    
    encoded_fps = int(video.get(cv2.CAP_PROP_FPS))
    real_fps = 1000  # The actual filming frame rate
    width = int(video.get(cv2.CAP_PROP_FRAME_WIDTH))
    height = int(video.get(cv2.CAP_PROP_FRAME_HEIGHT))
    total_frames = int(video.get(cv2.CAP_PROP_FRAME_COUNT))
    
    print(f"Encoded Video FPS: {encoded_fps}")
    print(f"Real FPS: {real_fps}")
    print(f"Video dimensions: {width}x{height}")
    print(f"Total frames: {total_frames}")
    print(f"Total data points: {len(data)}")

    original_data_hz = 100  # Original sensor data rate
    interpolated_data_hz = real_fps  # Interpolate to match video frame rate
    
    # Interpolate the pressure and acceleration magnitude data
    interpolated_pressure = interpolate_data(data['pres'], original_data_hz, interpolated_data_hz)
    interpolated_accmag = interpolate_data(data['accmag'], original_data_hz, interpolated_data_hz)
    
    # Calculate the offset between video frames and data points
    data_start_index = nadir_index * (interpolated_data_hz // original_data_hz) - nadir_frame - data_offset * (interpolated_data_hz // original_data_hz)
    
    fourcc = cv2.VideoWriter_fourcc(*'mp4v')
    out = cv2.VideoWriter(output_path, fourcc, encoded_fps, (width, height))
    
    for frame_count in range(total_frames):
        ret, frame = video.read()
        if not ret:
            break
        
        # Adjust the frame height
        frame = adjust_frame(frame, video_height_adjust)
        
        # Calculate the corresponding data index for the current frame
        data_index = data_start_index + frame_count
        
        if 0 <= data_index < len(interpolated_pressure):
            graph = create_graphs(interpolated_pressure, interpolated_accmag, data_index, interpolated_data_hz)
            frame = overlay_graph(frame, graph)
            
            # Add custom multi-line text to the top left of the video
            for i, line in enumerate(video_text):
                cv2.putText(frame, line, (10, 40 + i*40), 
                            cv2.FONT_HERSHEY_SIMPLEX, 1.0, (255, 255, 255), 2, cv2.LINE_AA)
        
        out.write(frame)
    
    video.release()
    out.release()
    
    print(f"Total frames processed: {total_frames}")

# Usage
video_path = './development/B69-0717094652_pres_vid_2024-07-17_09-46-45.mp4'
csv_path = './development/B69-0717094652_imp.csv'
output_path = './development/B69-0717094652_imp.mp4'
nadir_frame = 620  # Frame where pressure nadir occurs
nadir_index = 1443  # Row in CSV file where pressure nadir occurs
data_offset = 2    # Offset to align data correctly
video_height_adjust = 90  # Number of pixels to move the video up
video_text = ["PumpFlow & University of Hull",
    "TalTech RAPID 100hz",
    "Video speed: 1000fps",
    "Pump speed: 400 RPM"]
process_video(video_path, csv_path, output_path, nadir_frame, nadir_index, data_offset, video_height_adjust, video_text=video_text)
