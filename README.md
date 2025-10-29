# OS-Home-work
>[!NOTE]
> This project is developed under the **01204332 Operating System** course of **Department of Computer Engineering**, **Faculity of Engineering**, **Kasetsart University**.

> Developed by Nakharet Mueangphakdee (6610505438)

## Directory Structure
```bash
├── 1_parallel_6610505438
│   ├── parallel.py
│   └── requirements.txt
├── 2_cow_6610505438
│   ├── cow.py
│   └── requirements.txt
├── .gitignore
├── README.md
└── report_6610505438.pdf

```
## Installation prerequisite
### 1. OpenMPI
Linux 
```bash
sudo apt update
sudo apt install -y openmpi-bin libopenmpi-dev
```
After installation, verify:
```bash
mpirun --version
```

Windows

Download from official microsoft [OpenMPI](https://learn.microsoft.com/en-us/message-passing-interface/microsoft-mpi).
After installation, verify:
```bash
mpiexec --version
```
### 2. Python and pip
Linux 
```bash
sudo apt update
sudo apt install -y python3 pip
```
After installation, verify:
```bash
python --version
pip --version
```

Windows
Download from official python site [Python](https://www.python.org/downloads/)

After installation, verify:
```bash
python --version
pip --version
```
## How to run
### 1. Clone the Repository

```bash
git clone https://github.com/NareMu/os_homework_6610505438.git
cd ./os_homework_6610505438
```

### 2. Install Requirements
For parallel homework run 
```bash
pip install -r ./1_parallel_6610505438/requirements.txt
```

For cow homework run
```bash
pip install -r ./2_cow_6610505438/requirements.txt
```
### 3. Run the code
#### A) For parallel homework run by use open mpi

In Linux
```bash
mpirun -np <number_of_process> python ./1_parallel_6610505438/parallel.py <number_to_factorize>
```

In Windows
```bash
mpiexec -np <number_of_process> python .\1_parallel_6610505438\parallel.py <number_to_factorize>
```

Example
```bash
# In Linux example
mpirun -np 4 python ./1_parallel_6610505438/parallel.py 350243405507562291174415825999

# In Windows example
mpiexec -np 4 python .\1_parallel_6610505438\parallel.py 350243405507562291174415825999
```

#### B) For cow homework run
```bash
python ./2_cow_6610505438/cow.py <allocate size_in_MB> <modify_size_ratio(0.0-1.0)>
```
Example
```bash
python ./2_cow_6610505438/cow.py 120 0.8
```
