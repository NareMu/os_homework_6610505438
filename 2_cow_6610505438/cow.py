import numpy as np
import time
import os
import sys

def get_page_size() -> int:
    try:
        return os.sysconf('SC_PAGE_SIZE') 
    except:
        return 4096

def smaps_rollup(pid: int) -> dict[str,int]:
    path = f"/proc/{pid}/smaps_rollup"
    out = {}
    try:
        with open(path, "r", encoding="utf-8", errors="ignore") as f:
            for line in f:
                if ":" in line and "kB" in line:
                    key, val = line.split(":", 1)
                    num = val.strip().split()[0]
                    out[key.strip()] = int(num) * 1024  # to bytes
    except FileNotFoundError as e:
        print(f"Error: {e}")
    return out 

def cow_simulate(alloc_size_mb: int, modify_size_ratio: float ) -> None:
    page_size = get_page_size() 
    n_byte = alloc_size_mb * 2**20
    arr = np.random.randint(0, 256, size = n_byte, dtype = np.uint8)

    pid = os.fork()

    if pid == 0:
        child_pid = os.getpid()
        parent_pid = os.getppid()

        print("RSS before child edit the data")
        parent_info = smaps_rollup(parent_pid)
        child_info = smaps_rollup(child_pid)
        print(f"[Child] child process's RSS = {child_info["Rss"]/2**20} MiB parent process's RSS = {parent_info["Rss"]/2**20} MiB")
        print(f"[Child] child process's shared = {(child_info["Shared_Clean"]+child_info["Shared_Dirty"])/2**20} MiB parent process's shared = {(parent_info["Shared_Clean"]+parent_info["Shared_Dirty"])/2**20} MiB")
    
        time.sleep(0.5)
        n_modify = int(n_byte * modify_size_ratio)

        for i in range(0, n_modify, page_size):
            arr[i] = (int(arr[i]) + 1) % 256
        

        time.sleep(0.5)
        print("RSS after child edit the data")
        parent_info = smaps_rollup(parent_pid)
        child_info = smaps_rollup(child_pid)
        print(f"[Child] child process's RSS = {child_info["Rss"]/2**20} MiB parent process's RSS = {parent_info["Rss"]/2**20} MiB")
        print(f"[Child] child process's shared = {(child_info["Shared_Clean"]+child_info["Shared_Dirty"])/2**20} MiB parent process's shared = {(parent_info["Shared_Clean"]+parent_info["Shared_Dirty"])/2**20} MiB")
     
        print("[Child] child terminate")
        os._exit(0)
    else:
        os.waitpid(pid, 0)
        print("[Parent] child already terminate")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(f"Usage: python cow.py <alloc_size_mb> <modify_size_ratio>")
        sys.exit(1)
    
    try:
            # sys.argv stores arguments as strings, so we must convert them
            alloc_size = int(sys.argv[1])
            modify_ratio = float(sys.argv[2])

            # Validate ratio
            if not 0.0 <= modify_ratio <= 1.0:
                print(f"Error: modify_size_ratio must be between 0.0 and 1.0. Got: {modify_ratio}")
                sys.exit(1)
                
            # Call the main function
            cow_simulate(alloc_size_mb=alloc_size, modify_size_ratio=modify_ratio)

    except ValueError:
        print("Error: Invalid argument type.")
        print(f"  <alloc_size_mb> must be an integer ")
        print(f"  <modify_size_ratio> must be a float")
        sys.exit(1)

