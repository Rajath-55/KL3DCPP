import os
import subprocess
from math import ceil, log2


def get_number_of_nodes(filename):
  """
  Reads the first line of the file and returns the number of nodes.
  - If the first line is a power of 2, returns it directly.
  - Otherwise, returns the next power of 2.
  """
  with open(filename, 'r') as f:
    first_line = int(f.readline().strip())
  return first_line if (first_line & (first_line - 1) == 0) else 2 ** (ceil(log2(first_line)))


def run_command(command, filename, x):
  """
  Runs the command with x substituted and captures the cost value.
  """
  # Use f-strings for cleaner formatting (Python 3.6+)
  process = subprocess.Popen(f"{command}",
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
  output, _ = process.communicate()
  for line in output.decode().splitlines():
    if line.startswith("cost in CC"):
      return float(line.split()[-1])
  return None


def main():
  cpp_dict = {}
  python_results = {}
  for filename in os.listdir():
    if filename.startswith("Graph") and filename.endswith(".txt"):
      number_of_nodes = get_number_of_nodes(filename)
      for x in range(2, number_of_nodes + 1, 2):
        cpp_cost = run_command(f"./a.out {filename} 10 0 {x} | grep \"cost\"", filename, x)
        python_cost = run_command(f"python3 ./python_code/main_cmd.py {filename} 100 0 {x}", filename, x)
        if cpp_cost is not None:
          cpp_dict[f"layers={x}"] = cpp_cost
        if python_cost is not None:
          python_results[f"layers={x}"] = python_cost
  print("C++ Results:")
  print(cpp_dict)
  print("Python Results:")
  print(python_results)


if __name__ == "__main__":
  main()