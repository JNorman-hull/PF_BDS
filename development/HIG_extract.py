##two different methods for extracting raw binary values from HIG files


import struct

def examine_hig_file(filename, num_packets=5):
    with open(filename, 'rb') as file:
        for i in range(num_packets):
            packet = file.read(11)  # Read one packet (11 bytes)
            if not packet:
                break  # End of file

            # Unpack the data
            timestamp = struct.unpack('>i', packet[:4])[0]
            acc_x = struct.unpack('>h', packet[4:6])[0]
            acc_y = struct.unpack('>h', packet[6:8])[0]
            acc_z = struct.unpack('>h', packet[8:10])[0]
            end_byte = packet[10]

            # Print the data
            print(f"Packet {i+1}:")
            print(f"  Timestamp: {timestamp}")
            print(f"  Acc X: {acc_x}")
            print(f"  Acc Y: {acc_y}")
            print(f"  Acc Z: {acc_z}")
            print(f"  End byte: {end_byte:02x}")
            print(f"  Raw bytes: {packet.hex()}")
            print()

# Usage
examine_hig_file('./RAW_data/RAPID/B11-0529102932.HIG')


def print_raw_binary_data(filename: str, num_packets: int = 5):
    packetSize = 11  # Define the packet size for .HIG files
    with open(filename, 'rb') as file_ID:
        for i in range(num_packets):
            packet = file_ID.read(packetSize)
            if len(packet) < packetSize:
                break  # Stop if we reach the end of the file or packet is incomplete
            time_raw = struct.unpack('>I', packet[:4])[0]
            acc_x = struct.unpack('>h', packet[4:6])[0]
            acc_y = struct.unpack('>h', packet[6:8])[0]
            acc_z = struct.unpack('>h', packet[8:10])[0]
            print(f"Packet {i+1}: Time: {time_raw}, AccX: {acc_x}, AccY: {acc_y}, AccZ: {acc_z}")

# Usage
print_raw_binary_data('./RAW_data/RAPID/B11-0529102932.HIG', num_packets=10)
