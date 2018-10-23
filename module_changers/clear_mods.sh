# .bashrc

# Source global definitions
if [ -f /etc/bashrc ]; then
	. /etc/bashrc
fi

# Script to remove all loaded modules and direct output to stderr
# Would not clear mods properly so made this file to run on the cluster...sorry

echo "Clearing modules" >&2
module purge >&2
echo "modules cleared" >&2
module list
