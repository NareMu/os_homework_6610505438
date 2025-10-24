# OS-Home-work
> This project is developed under the **01204332 Operating System** course of **Department of Computer Engineering**, **Faculity of Engineering**, **Kasetsart University**.

> This project for **educational purpose** only.
## Directory Structure
```hs
os_homework_6610505438
|-- 1_parallel_6610505438
|   └--parallel.py
|-- README.md
└-- report_6610505438.pdf
```
## Installation prerequisite
### 1. Install OpenMPI
Linux (Ubuntu / Debian)
```bash
sudo apt update
sudo apt install -y openmpi-bin libopenmpi-dev
```
After installation, verify:
```bash
mpirun --version
```
macOS (using Homebrew)

```bash
brew install open-mpi
```

After installation, verify:
```bash
mpirun --version
```
Windows

Download from official microsoft [OpenMPI](https://learn.microsoft.com/en-us/message-passing-interface/microsoft-mpi).


## How to run
### 1. Clone the Repository

```bash
git clone https://github.com/NareMu/os_homework_6610505438.git
cd ./os_homework_6610505438
```

### 2. Install Requirements
```bash
pip install -r requirements.txt
```
### 3. Run the code
For parallel homework run by use open mpi

In Mac/Linux
```bash
mpirun -np {number of process} python ./1_parallel_6610505438/parallel.py {number to factorize}
```

In Window
```bash
mpiexec -np {number of process} python ./1_parallel_6610505438/parallel.py {number to factorize}
```
