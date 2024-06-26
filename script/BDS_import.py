#Python BDS txt > CSV conversion
#Adapated from https://github.com/jtuhtan/RAPID/tree/main

import struct
from abc import ABC, abstractmethod
from pathlib import Path
from jupyter_client import BlockingKernelClient
from typing import Tuple
from datetime import datetime
import quaternion

from matplotlib import pyplot as plt
plt.style.use("seaborn-v0_8-whitegrid")

import numpy as np
import pandas as pd

class Rapid(ABC):
    def __init__(self, filename: str) -> None:  # Attributes same for all Sensors can be placed here
        """This is an abstract class (ABC).
        It provides methods which are the same between the sensor / subclasses.
        This class is never called by the user, only by its children.

        """
        self.filename = Path(filename)
        self.dir_csv = Path("csv/")
        self.dir_plots = Path("plots/")

    def _mkdir(self, path_object: Path) -> None:
        path_object.mkdir(parents=True, exist_ok=True)

    def _save_as_csv(self, data, **kwargs) -> None:
        self._mkdir(self.dir_csv)
        data.to_csv(
            (self.dir_csv / self.filename.stem).with_suffix(".csv"),
            sep=",",
            index=False,
            **kwargs
        )

    def _process_and_save(self, savecsv, **kwargs):
        data = self._read_data()
        data, summary_info = self._post_process(data)
        if savecsv:
            self._save_as_csv(data, **kwargs)
        return data, summary_info

    def _read_data(self) -> pd.DataFrame:
        with open(self.filename.as_posix(), mode="r+b") as f:
            binary_data = f.read()
        len_fmt = struct.calcsize(self.fmt)
        max_accept_len = (len(binary_data) // len_fmt) * len_fmt
        bin_new = binary_data[0:max_accept_len]
        iter = struct.iter_unpack(self.fmt, bin_new)
        data = [x for x in iter]
        data = pd.DataFrame(data, columns=self.column_names_raw)
        return data
      
    def parse_filename_info(self) -> dict:
        """Extracts the sensor name, date, and time from the filename."""
        filename = self.filename.stem  # Use stem to avoid the file extension

        # Assuming a filename format like 'C010706133321.txt'
        sensor = filename[:3]  # First three characters are the sensor name
        date_str = filename[3:7]  # Characters 4 to 7 are the date in MMDD format
        time_str = filename[7:13]  # Characters 8 to 13 are the time in HHMMSS format

        date_deploy = datetime.strptime(date_str, "%m%d").strftime("%d/%m")
        time_deploy = f"{time_str[:2]}:{time_str[2:4]}:{time_str[4:]}"

        return {
            'sensor': sensor,
            'date_deploy': date_deploy,
            'time_deploy': time_deploy
        }


    @abstractmethod
    def _class_specific_config(self) -> Tuple[int, list[str]]:
        """Provide summary interval and relevant fields."""
        pass

    def check_unchanged_sensors(self, data: pd.DataFrame, threshold: int = 96) -> dict[str, bool]:
        """Check if any sensor data remains constant over a given threshold.
    
        Parameters:
        - data (pd.DataFrame): The dataset containing sensor readings.
        - threshold (int): The number of rows to check for unchanged data.
    
        Returns:
        - list[str]: List of sensors that have remained unchanged.
        """
        #threhold @ 96 = check every 1s
        #P1 = Left, P2 = Center, P3 = Right
        unchanged_sensors = {}
        for sensor in ["P1", "P2", "P3"]:
            if sensor in data.columns:
                series = data[sensor]
                unchanged = series.rolling(threshold).apply(lambda x: pd.Series(x).nunique() == 1, raw=True).max() == 1
                unchanged_sensors[sensor] = not unchanged
        return unchanged_sensors
      
    def compute_pres(self, data: pd.DataFrame, valid_sensors: dict[str, bool]) -> pd.DataFrame:
        """Compute the 'pres' value based on the working sensors.
    
        Parameters:
        - data (pd.DataFrame): The dataset containing sensor readings.
        - valid_sensors (Dict[str, bool]): Dictionary indicating which sensors are working.
    
        Returns:
        - pd.Series: Computed 'pres' values.
        """
        # Define fallback logic for computing the 'pres' value
        #If only one sensor works, return the value, if more than one works, return average
        def compute_row(row):
            active_sensors = [sensor for sensor, is_valid in valid_sensors.items() if is_valid]
            values = [row[sensor] for sensor in active_sensors]
            if len(values) == 0:
                return [np.nan, ','.join(active_sensors)]
            elif len(values) == 1:
                return [values[0], ','.join(active_sensors)]
            else:
                return [np.mean(values), ','.join(active_sensors)]
    
        result = data.apply(compute_row, axis=1, result_type='expand')
        return pd.DataFrame(result.values, columns=['pres', 'sensors_used'])  

    def count_threshold_exceedances(self, data: pd.DataFrame, column: str, threshold: float) -> int:
        """Count exceedances of a specified threshold in the given column."""
        filtered = data[(data[column] >= threshold)]
        return filtered.shape[0]

    def _post_process(self, data: pd.DataFrame) -> Tuple[pd.DataFrame, dict[str, float]]:
        """Common post-processing logic shared between subclasses."""
        data["time"] = (data["time"] - data["time"].iloc[0]) / 1000
        
        if data["time"].max() >= 5000:
            time_warning = "ERROR: Time series incorrect. Data cannot be used."
        else:
            time_warning = ""
            
        valid_sensors = self.check_unchanged_sensors(data)
        pres_data = self.compute_pres(data, valid_sensors)
        sample_speed, relevant_fields = self._class_specific_config()
        
        data.insert(1, "pres", pres_data['pres'])
        data.insert(1, "accmag", np.linalg.norm(data[["accx", "accy", "accz"]], axis=-1)) #computes the Euclidean norm (magnitude) across the acceleration components along the x, y, and z axes.
        data["accmag"] -= 9.81 #removes earth gravity
        data = data[relevant_fields]

        #Calculate duration of sensor deployment
        #where sample_speed = 96 rows for 1s
        num_seconds = len(data)// sample_speed
        minutes, seconds = divmod(num_seconds, 60)
        duration = f"{minutes:02}:{seconds:02}"

        #consider 100 m/s2 as an impact indicator
        max_acc_g_force = data['accmag'].max()
        if max_acc_g_force >= 100:
            acc_warning = "ACC: Acceleration event >= 100 m/s2 detected"
        else:
            acc_warning = ""

        pres_min_index = data['pres'].idxmin()
        pres_min_time = data['time'][pres_min_index]
        pres_max_index = data['pres'].idxmax()
        pres_max_time = data['time'][pres_max_index]
        acc_max_index = data['accmag'].idxmax()
        acc_max_time = data['time'][acc_max_index]
        
        unchanged_sensors = [sensor for sensor, is_valid in valid_sensors.items() if not is_valid]
        sensors_used_summary = pres_data['sensors_used'].unique()
        sensors_used_summary_text = ', '.join(sensors_used_summary)

        if unchanged_sensors:
            warning_message = f"WARNING: Pressure data error ({', '.join(unchanged_sensors)}). Average taken from {sensors_used_summary_text}"
        else:
            warning_message = "All pressure sensors working"

        if time_warning:
            warning_message = time_warning + "; " + warning_message
            
        if acc_warning:
            warning_message = acc_warning + "; " + warning_message

        summary_info = {
            'duration[mm:ss]': duration,
            'pres_min[mbar]': data['pres'].min(),
            'pres_min[time]': pres_min_time,
            'pres_max[mbar]': data['pres'].max(),
            'pres_max[time]': pres_max_time,
            'acc_max[m/s2]': data['accmag'].max(),
            'acc_max[time]': acc_max_time,
            'messages': warning_message
        }

        return data, summary_info
      
    def plot_data_overview(self, save: bool = True, show: bool = False) -> None:
        """Plots an overview for the generated data.
        This is primarily to spot problems before further user-processing.
        The save-option can be helpful as a visual aid for the future as to
        which measurements measured what.

        Parameters
        ----------
        save : bool
            Save the file at location specified in <sensorclass>.dir_plots, by default True
        show : bool
            Show an interactive plot when executed, by default False
        """
        # Determine if the data is from IMP format
        is_imp = "Accel_Mag (g)" in self.data.columns
        
        t = self.data["time"][::10]
        
        if is_imp:
            pres = self.data["Pressure (mbar)"].rolling(10).mean()[::10]
            accmag = self.data["Accel_Mag (g)"].rolling(10).mean()[::10]
        else:
            pres = self.data["pres"].rolling(10).mean()[::10]
            accmag = self.data["accmag"].rolling(10).mean()[::10]
        
        color = "C0"
        fig, ax1 = plt.subplots(figsize=(25, 5))
        ax1.set_xlabel("time [s]")
        ax1.set_ylabel("Pressure [hPa]", color=color)
        ax1.plot(t, pres, color=color)
        ax1.tick_params(axis="y", labelcolor=color)
        ax1.ticklabel_format(useOffset=False)

        ax2 = ax1.twinx()
        color = "C1"
        ax2.set_ylabel("Acceleration Magnitude (m/s\u00b2)", color=color)
        ax2.plot(t, accmag, color=color)
        ax2.tick_params(axis="y", labelcolor=color)
        fig.tight_layout()

        if save == True:
            self._mkdir(self.dir_plots)
            new_filename =  self.filename.stem + "_pres_acc"
            plt.savefig((self.dir_plots / new_filename).with_suffix(".png"))
        if show == True:
            plt.show()
        plt.close()
        
    def plot_pressure_overview(self, save: bool = True, show: bool = False) -> None:
        """Plot individual raw values of P1, P2, and P3 to identify malfunctioning sensors."""
        t = self.data["time"][::10]  # Downsampling for visualization
        p1 = self.data["P1"].rolling(10).mean()[::10]  # Rolling average for smoothing
        p2 = self.data["P2"].rolling(10).mean()[::10]
        p3 = self.data["P3"].rolling(10).mean()[::10]
    
        fig, ax = plt.subplots(figsize=(25, 5))
        ax.set_xlabel("time [s]")
        ax.set_ylabel("Pressure [hPa]")
    
        # Plot individual pressure sensors
        ax.plot(t, p1, color="C1", label="P1") #left sensor
        ax.plot(t, p2, color="C2", label="P2") #center sensor
        ax.plot(t, p3, color="C3", label="P3") #right sensor
        ax.legend()
        
        #check for unchanged pressure values and report error
        valid_sensors = self.check_unchanged_sensors(self.data)
        pres_data = self.compute_pres(self.data, valid_sensors)
        unchanged_sensors = [sensor for sensor, is_valid in valid_sensors.items() if not is_valid]
        sensors_used_summary = pres_data['sensors_used'].unique()
        sensors_used_summary_text = ', '.join(sensors_used_summary)
        if unchanged_sensors:
          warning_message = f"WARNING: Pressure data error ({', '.join(unchanged_sensors)}). Average taken from {sensors_used_summary_text}"
        else:
          warning_message = "All pressure sensors working"
        
        if warning_message:
          ax.text(0.5, 0.9, warning_message, fontsize=25, color="red", ha="center", transform=ax.transAxes)
        
        fig.tight_layout()

        if save == True:
            new_filename = self.filename.stem + "_pres_validation" 
            plt.savefig((self.dir_plots / new_filename).with_suffix(".png"))
        if show == True:
            plt.show()
        plt.close() 
        

class BDS100(Rapid):
    def __init__(self, filename: str, savecsv: bool = True, **kwargs) -> None:
        """This class processes BDS measurements at 100 Hz. 
        At that speed, the internal data fusion algorithm computes an absolute orientation, 
        which is used to get from a relative to an absolute acceleration.
        Acceleration magnitude is included in the ouput and 
        has gravity already removed.

        Parameters
        ----------
        filename : str
            Relative file location + filename of the measurement
        savecsv : bool, optional
            Saves the processed data as a csv file if True, by default True
        **kwargs : optional
            Keyword arguments for changing how the csv is generated and are feeded directly into pd.read_csv().
        """
        super().__init__(filename)
        self.dir_csv = self.dir_csv / "BDS100"
        self.dir_plots = self.dir_plots / "BDS100"
        self.fmt = "HI22f4B"  # format string to set byteorder
        self.class_name = "100hz"
        self.column_names_raw = [
            "sample rate", "time", "P1", "T1", "P2", "T2", "P3", "T3", "eul head", "eul roll", "eul pitch", "quat w",
            "quatx", "quaty", "quatz", "magx", "magy", "magz", "accx", "accy", "accz", "gyrox", "gyroy", "gyroz", "calmag",
            "calacc", "calgyro", "calimu"]
        self.data, _ = super()._process_and_save(savecsv, **kwargs)

    def _class_specific_config(self) -> Tuple[int, list[str]]:
        return 96, ["time", "pres", "P1", "P2", "P3", "accmag", "magx", "magy", "magz", "quat w", "quatx", "quaty", "quatz", "accx", "accy", "accz"]
  #96 rows = 1s
  
    def _post_process(self, data: pd.DataFrame) -> Tuple[pd.DataFrame, dict[str, float]]:
        """Override post_process to add absolute acceleration calculations."""
        # Perform the base class processing
        data, summary_info = super()._post_process(data)
        # Add absolute acceleration calculations
        data = self._absolute_orientation(data)
        relevant_fields = ["time", "pres", "P1", "P2", "P3", "accmag", "magx", "magy", "magz", "absaccx", "absaccy", "absaccz"]
        data = data[relevant_fields]
        return data, summary_info
   
    def _absolute_orientation(self, data: pd.DataFrame) -> pd.DataFrame:
        # Translate body acc with earths mag field to abs reference frame
        quat_ref_frame = quaternion.as_quat_array(
            data[["quat w", "quatx", "quaty", "quatz"]]
        )

        acc_earth = np.zeros((0, 3))
        for idx, q in enumerate(quat_ref_frame):
            acc = quaternion.rotate_vectors(q, data.loc[idx, "accx":"accz"])
            acc_earth = np.vstack((acc_earth, acc))
        acc_earth[:, 2] -= 9.81
        data.insert(5, "absaccx", acc_earth[:, 0])
        data.insert(6, "absaccy", acc_earth[:, 1])
        data.insert(7, "absaccz", acc_earth[:, 2])
        
        return data
      
class BDS250(Rapid):
    def __init__(self, filename: str, savecsv: bool = True, **kwargs) -> None:
        """This class processes BDS measurements at 250 Hz. 
        At that speed, there is no absolute orientation computation, 
        outputs are absolute pressure, accelerometer and gyroscope.
        Acceleration magnitude is included in the ouput and 
        has gravity already removed.

        Parameters
        ----------
        filename : str
            Relative file location + filename of the measurement
        savecsv : bool, optional
            Saves the processed data as a csv file if True, by default True
        **kwargs : optional
            Keyword arguments for changing how the csv is generated and are feeded directly into pd.read_csv().
        """
        super().__init__(filename)
        self.dir_csv = self.dir_csv / "BDS250"
        self.dir_plots = self.dir_plots / "BDS250"
        self.fmt = "HI12f4B"  # format string to set byteorder
        self.class_name = "250hz"
        self.column_names_raw = [
            "samplerate", "time", "P1", "T1", "P2", "T2", "P3", "T3", "accx", "accy", "accz", "gyrox", "gyroy",
            "gyroz", "calmag", "calacc", "calgyro", "calimu"]
        self.data, _ = super()._process_and_save(savecsv, **kwargs)

    def _class_specific_config(self) -> Tuple[int, list[str]]:
        return 96, ["time", "pres", "P1", "P2", "P3", "accmag"]
    #96 rows = 1s. Even when sensor set to 250hz, data log is at 100hz
    
class IMP(Rapid):
    def __init__(self, filename: str, savecsv: bool = True, **kwargs) -> None:
        super().__init__(filename)
        self.dir_csv = self.dir_csv / "RAPID_IMP"
        self.dir_plots = self.dir_plots / "RAPID_IMP"
        self.class_name = "100_imp"
        self.packetSize = 29
        self.FS = 2000
        self.IMU_PREC = 3
        self.P_PREC = 1
        self.T_BAT_PREC = 2
        self.gain_ac = 0.005 #/ 9.81
        self.gain_gy = 0.1
        self.gain_mg = 0.1
        self.gain_pr = 0.1
        self.gain_t = 0.01
        self.gain_bt = 0.01
        self.column_names_raw = [
            'Time (s)', 'Accel_X (g)', 'Accel_Y (g)', 'Accel_Z (g)', 'Accel_Mag (g)',
            'Gyro_X (deg/s)', 'Gyro_Y (deg/s)', 'Gyro_Z (deg/s)', 'Mag_X (mT)',
            'Mag_Y (mT)', 'Mag_Z (mT)', 'Pressure (mbar)', 'P_Temp (C)', 'Battery (V)'
        ]
        self.data, _ = self._process_and_save(savecsv, **kwargs)

    def _class_specific_config(self) -> Tuple[int, list[str]]:
        return 96, self.column_names_raw

    def _read_data(self) -> pd.DataFrame:
        with open(self.filename.as_posix(), 'rb') as file_ID:
            fstat = os.stat(self.filename)
            flen = (fstat.st_size // self.packetSize) - 1

            TimeRaw = []
            TimeSpot = []

            for _ in range(flen):
                time_raw = struct.unpack('>i', file_ID.read(4))[0]
                TimeRaw.append(time_raw)
                TimeSpot.append(file_ID.tell())
                file_ID.seek(25, 1)

            DataRaw = np.zeros((flen, 12), dtype=np.int16)
            DataRawP = np.zeros(flen, dtype=np.uint16)

            for it in range(flen):
                if it == 0:
                    file_ID.seek(4, 0)
                else:
                    DataRaw[it, 0] = struct.unpack('>h', file_ID.read(2))[0]
                    DataRaw[it, 1] = struct.unpack('>h', file_ID.read(2))[0]
                    DataRaw[it, 2] = struct.unpack('>h', file_ID.read(2))[0]
                    DataRaw[it, 3] = struct.unpack('>h', file_ID.read(2))[0]
                    DataRaw[it, 4] = struct.unpack('>h', file_ID.read(2))[0]
                    DataRaw[it, 5] = struct.unpack('>h', file_ID.read(2))[0]
                    DataRaw[it, 6] = struct.unpack('>h', file_ID.read(2))[0]
                    DataRaw[it, 7] = struct.unpack('>h', file_ID.read(2))[0]
                    DataRaw[it, 8] = struct.unpack('>h', file_ID.read(2))[0]
                    file_ID.seek(2, 1)
                    DataRaw[it, 9] = struct.unpack('>h', file_ID.read(2))[0]
                    DataRaw[it, 10] = struct.unpack('>h', file_ID.read(2))[0]
                    file_ID.seek(5, 1)

            DataRaw[0, :] = DataRaw[1, :]

            for it in range(flen):
                if it == 0:
                    file_ID.seek(TimeSpot[0] + 4 + (2 * 7), 0)
                else:
                    DataRawP[it] = struct.unpack('>H', file_ID.read(2))[0]
                    file_ID.seek(27, 1)

            DataRawP[0] = DataRawP[1]

            RAPIDIMP = {
                'td': np.array(TimeRaw, dtype=np.float64),
                'ts': np.array(TimeRaw, dtype=np.float64) / self.FS,
                'ax': np.round(DataRaw[:, 0] * self.gain_ac, self.IMU_PREC),
                'ay': np.round(DataRaw[:, 1] * self.gain_ac, self.IMU_PREC),
                'az': np.round(DataRaw[:, 2] * self.gain_ac, self.IMU_PREC),
                'gx': np.round(DataRaw[:, 3] * self.gain_gy, self.IMU_PREC),
                'gy': np.round(DataRaw[:, 4] * self.gain_gy, self.IMU_PREC),
                'gz': np.round(DataRaw[:, 5] * self.gain_gy, self.IMU_PREC),
                'mx': np.round(DataRaw[:, 6] * self.gain_mg, self.IMU_PREC),
                'my': np.round(DataRaw[:, 7] * self.gain_mg, self.IMU_PREC),
                'mz': np.round(DataRaw[:, 8] * self.gain_mg, self.IMU_PREC),
                'p': np.round(DataRawP * self.gain_pr, self.P_PREC),
                't': np.round(DataRaw[:, 9] * self.gain_t, self.T_BAT_PREC),
                'b': np.round(DataRaw[:, 10] * self.gain_bt, self.T_BAT_PREC)
            }

            aMag = np.round(np.sqrt(RAPIDIMP['ax'] ** 2 + RAPIDIMP['ay'] ** 2 + RAPIDIMP['az'] ** 2) -9.81, self.IMU_PREC)
            
            dataExportCSV = np.column_stack((
                RAPIDIMP['ts'], RAPIDIMP['ax'], RAPIDIMP['ay'], RAPIDIMP['az'], aMag,
                RAPIDIMP['gx'], RAPIDIMP['gy'], RAPIDIMP['gz'], RAPIDIMP['mx'], RAPIDIMP['my'],
                RAPIDIMP['mz'], RAPIDIMP['p'], RAPIDIMP['t'], RAPIDIMP['b']
            ))

            return pd.DataFrame(dataExportCSV, columns=self.column_names_raw)

    def _post_process(self, data: pd.DataFrame) -> Tuple[pd.DataFrame, dict[str, float]]:
        # Calculate duration
        data["time"] = data["Time (s)"]
        num_seconds = len(data) // 100
        minutes, seconds = divmod(num_seconds, 60)
        duration = f"{int(minutes):02}:{int(seconds):02}"

        # Find min and max pressure times
        pres_min_index = data['Pressure (mbar)'].idxmin()
        pres_min_time = data['Time (s)'][pres_min_index]
        pres_max_index = data['Pressure (mbar)'].idxmax()
        pres_max_time = data['Time (s)'][pres_max_index]

        # Find max acceleration time
        acc_max_index = data['Accel_Mag (g)'].idxmax()
        acc_max_time = data['Time (s)'][acc_max_index]
        
        if data["time"].max() >= 5000:
            time_warning = "ERROR: Time series incorrect. Data cannot be used."
        else:
            time_warning = ""
            
        #consider 100 m/s2 as an impact indicator
        max_acc_g_force = data['Accel_Mag (g)'].max()
        if max_acc_g_force >= 100:
            acc_warning = "ACC: Acceleration event >= 100 m/s2 detected"
        else:
            acc_warning = ""

        warning_message = "No errors detected"
        if time_warning or acc_warning:
            warning_message = "; ".join(filter(None, [time_warning, acc_warning]))
            
        # Summary information specific to IMP
        summary_info = {
            'duration[mm:ss]': duration,
            'pres_min[mbar]': data['Pressure (mbar)'].min(),
            'pres_min[time]': pres_min_time,
            'pres_max[mbar]': data['Pressure (mbar)'].max(),
            'pres_max[time]': pres_max_time,
            'acc_max[m/s2]': data['Accel_Mag (g)'].max(),
            'acc_max[time]': acc_max_time,
            'messages': warning_message
        }
        return data, summary_info
      
    def parse_filename_info(self) -> dict:
        """Extracts the sensor name, date, and time from the filename."""
        filename = self.filename.stem  # Use stem to avoid the file extension

       # Assuming a filename format like 'B80-0612182322' for IMP files
        if '-' in filename:
            sensor, date_time = filename.split('-')
        else:
            sensor = filename[:3]
            date_time = filename[3:]

        date_str = date_time[:4]  # Characters 1 to 6 are the date in DDMMYY format
        time_str = date_time[4:]  # Characters 7 to 12 are the time in HHMMSS format

        date_deploy = datetime.strptime(date_str, "%m%d").strftime("%d/%m")
        time_deploy = f"{time_str[:2]}:{time_str[2:4]}:{time_str[4:]}"

        return {
            'sensor': sensor,
            'date_deploy': date_deploy,
            'time_deploy': time_deploy
        }
        
class HIG(Rapid):
    def __init__(self, filename: str, savecsv: bool = True, **kwargs) -> None:
        super().__init__(filename)
        self.dir_csv = self.dir_csv / "RAPID_HIG"
        self.dir_plots = self.dir_plots / "RAPID_HIG"
        self.class_name = "100_hig"
        self.packetSize = 11
        self.FS = 2000
        self.HIG_PREC = 1
        self.gain_hig = 0.1
        self.column_names_raw = [
            'Time (s)', 'HIGAccel_X (g)', 'HIGAccel_Y (g)', 'HIGAccel_Z (g)', 'HIGAccel_Mag (g)'
        ]
        self.data, _ = self._process_and_save(savecsv, **kwargs)

    def _class_specific_config(self) -> Tuple[int, list[str]]:
        return 2000, self.column_names_raw

    def _read_data(self) -> pd.DataFrame:
        with open(self.filename.as_posix(), 'rb') as file_ID:
            fstat = os.stat(self.filename)
            flen = fstat.st_size // self.packetSize

            TimeRaw = []
            DataRaw = np.zeros((flen, 3), dtype=np.int16)
            
            for it in range(flen):
                packet = file_ID.read(self.packetSize)
                
                time_raw = struct.unpack('>I', packet[:4])[0]
                TimeRaw.append(time_raw)

                DataRaw[it, 0] = struct.unpack('>h', packet[4:6])[0]  # acc X
                DataRaw[it, 1] = struct.unpack('>h', packet[6:8])[0]  # acc Y
                DataRaw[it, 2] = struct.unpack('>h', packet[8:10])[0]  # acc Z
                
            TimeRaw = np.array(TimeRaw, dtype=np.float64)
            ts = TimeRaw / self.FS
            ax = np.round(DataRaw[:, 0] * self.gain_hig, self.HIG_PREC)
            ay = np.round(DataRaw[:, 1] * self.gain_hig, self.HIG_PREC)
            az = np.round(DataRaw[:, 2] * self.gain_hig, self.HIG_PREC)
            aMag = np.round(np.sqrt(ax ** 2 + ay ** 2 + az ** 2), self.HIG_PREC)

            dataExportCSV = np.column_stack((ts, ax, ay, az, aMag))
            return pd.DataFrame(dataExportCSV, columns=self.column_names_raw)
          
    def _post_process(self, data: pd.DataFrame) -> Tuple[pd.DataFrame, dict[str, float]]:
        # Calculate duration
        data["time"] = data["Time (s)"]
        num_seconds = len(data) / self.FS
        minutes, seconds = divmod(num_seconds, 60)
        duration = f"{int(minutes):02}:{int(seconds):02}"

        # Find max acceleration time
        acc_max_index = data['HIGAccel_Mag (g)'].idxmax()
        acc_max_time = data['Time (s)'][acc_max_index]

        if data["time"].max() >= 5000:
            time_warning = "ERROR: Time series incorrect. Data cannot be used."
        else:
            time_warning = ""

        # Check for acceleration warning
        max_acc_g_force = data['HIGAccel_Mag (g)'].max()
        if max_acc_g_force >= 100:
            acc_warning = "HIG event >= 100 g detected"
        else:
            acc_warning = ""

        # Warning message
        warning_message = ""
        if time_warning or acc_warning:
            warning_message = "; ".join(filter(None, [time_warning, acc_warning]))

        # Summary information specific to HIG
        summary_info = {
            'duration[mm:ss]': duration,
            'HIG_max[g]': max_acc_g_force,
            'HIG_max[time]': acc_max_time,
            'messages': warning_message
        }
        return data, summary_info
      
    def parse_filename_info(self) -> dict:
        """Extracts the sensor name, date, and time from the filename."""
        filename = self.filename.stem  # Use stem to avoid the file extension

       # Assuming a filename format like 'B80-0612182322' for HIG files
        if '-' in filename:
            sensor, date_time = filename.split('-')
        else:
            sensor = filename[:3]
            date_time = filename[3:]

        date_str = date_time[:4]  # Characters 1 to 6 are the date in DDMMYY format
        time_str = date_time[4:]  # Characters 7 to 12 are the time in HHMMSS format

        date_deploy = datetime.strptime(date_str, "%m%d").strftime("%d/%m")
        time_deploy = f"{time_str[:2]}:{time_str[2:4]}:{time_str[4:]}"

        return {
            'sensor': sensor,
            'date_deploy': date_deploy,
            'time_deploy': time_deploy
        }

if __name__ == "__main__":
    current_date = datetime.now().strftime("%d%m%y")
    summary_data = []
    n_files_w_pres_errors = 0
    n_files_w_time_errors = 0
    n_files_w_acc = 0
    n_files_w_hig = 0

#Process BDS100 files
    bds100_files = list(Path("./RAW_data/BDS_100").glob("*.txt"))
    print(f"Processing {len(bds100_files)} BDS100 files...")
    for f in bds100_files:    
        filename = f
        print(f"Processing {filename.name}")
        mymeas = BDS100(filename)
        _, summary_info = mymeas._process_and_save(True)
        mymeas.plot_data_overview()
        # mymeas.plot_acc_mag_overview()
        mymeas.plot_pressure_overview()
        
        if "WARNING" in summary_info["messages"]:
            n_files_w_pres_errors += 1
          
        if "ERROR" in summary_info["messages"]:
            n_files_w_time_errors += 1
            
        if "ACC" in summary_info["messages"]:
            n_files_w_acc += 1
        
        file_info = mymeas.parse_filename_info()
        summary_data.append({
            'file': filename.stem,
            'class': mymeas.class_name,
            **file_info,
            **summary_info
        })
        print(f"{filename.name} complete")


#Process BDS250 files
    bds250_files = list(Path("./RAW_data/BDS_250").glob("*.txt"))
    print(f"Processing {len(bds250_files)} BDS250 files...")
    for f in bds250_files:
        filename = f
        print(f"Processing {filename.name}")
        mymeas = BDS250(filename)
        _, summary_info = mymeas._process_and_save(True)
        mymeas.plot_data_overview()
        mymeas.plot_pressure_overview()
        
        if "WARNING" in summary_info["messages"]:
            n_files_w_pres_errors += 1
        
        if "ERROR" in summary_info["messages"]:
            n_files_w_time_errors += 1
        
        if "ACC" in summary_info["messages"]:
            n_files_w_acc += 1
        
        file_info = mymeas.parse_filename_info()
        summary_data.append({
            'file': filename.stem,
            'class': mymeas.class_name,
            **file_info,
            **summary_info
        })
        print(f"{filename.name} complete")
        
# Process IMP files
    imp_files = list(Path("./RAW_data/RAPID").glob("*.IMP"))
    print(f"Processing {len(imp_files)} RAPID IMP files...")
    for f in imp_files:
        filename = f
        print(f"Processing {filename.name}")
        mymeas = IMP(filename)
        _, summary_info = mymeas._process_and_save(True)
        mymeas.plot_data_overview()

        if "ERROR" in summary_info["messages"]:
            n_files_w_time_errors += 1
        
        if "ACC" in summary_info["messages"]:
            n_files_w_acc += 1
    
        file_info = mymeas.parse_filename_info()
        summary_data.append({
            'file': filename.stem,
            'class': mymeas.class_name,
            **file_info,
            **summary_info
        })
        print(f"{filename.name} complete")
        
        # Process HIG files
    hig_files = list(Path("./RAW_data/RAPID").glob("*.HIG"))
    print(f"Processing {len(hig_files)} RAPID HIG files...")
    for f in hig_files:
        filename = f
        print(f"Processing {filename.name}")
        mymeas = HIG(filename)
        _, summary_info = mymeas._process_and_save(True)

        if "ERROR" in summary_info["messages"]:
            n_files_w_time_errors += 1
        
        if "HIG" in summary_info["messages"]:
            n_files_w_hig += 1

        file_info = mymeas.parse_filename_info()
        summary_data.append({
        'file': filename.stem,
        'class': mymeas.class_name,
        **file_info,
        **summary_info
        })
        print(f"{filename.name} complete")

#Create summary CSV
    summary_df = pd.DataFrame(summary_data)
    summary_csv_filename = f"{current_date}_batch_summary.csv"
    summary_csv_path = Path("csv") / summary_csv_filename
    summary_df.to_csv(summary_csv_path, index=False)
    
    print("Batch sensor processing complete.")
    print(f"{len(bds100_files) + len(bds250_files) + len(imp_files) + len(hig_files)} total sensor files processed")
    print(f"{len(bds100_files) + len(bds250_files)} BDS files processed")
    print(f"{len(imp_files) + len(hig_files)} RAPID sensor files processed")
    print(f"{n_files_w_pres_errors}/{len(bds100_files) + len(bds250_files)} BDS files contain pressure errors")
    print(f"{n_files_w_time_errors}/{len(summary_data)} sensor files contain time series errors")
    print(f"{n_files_w_acc}/{len(bds100_files) + len(bds250_files) + len(imp_files)} BDS/RAPID IMU contain acceleration >= 100 m/s2")
    print(f"{n_files_w_hig}/{len(hig_files)} RAPID HIG contain events >= 100g")
