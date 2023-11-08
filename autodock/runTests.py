# General imports
import subprocess, sys

def runTests():
	"""
	This function runs the runTests.py script inside package pwchem.
	Note: To run this script, the scipion3 env must be active.
	"""
	# Define the command to run the script
	command = ["python", "-m", "pwchem.runTests"] + sys.argv[1:]

	# Run the script as a separate process
	try:
		subprocess.run(command, check=True)
	except subprocess.CalledProcessError:
		sys.exit(1)

# Call to main execution
if __name__ == "__main__":
	runTests()
