
# Functions to log to the screen and to a file
# All of them expect a file as the first argument and a message as the second
#  10/01/2020: Fabio Almeida. 

log() {
	echo -e "[INFO] \033[1;34m${FUNCNAME[1]}:\033[m $2"
	echo "$(date -u '+%d-%m-%Y %T %Z') [INFO] ${FUNCNAME[1]}: $2" >> "$1"
}

log_warn() {
	echo -e "[\033[1;33mWARN\033[m] \033[1;34m${FUNCNAME[1]}:\033[m $2" >&2
	echo "$(date -u '+%d-%m-%Y %T %Z') [WARN] ${FUNCNAME[1]}: $2" >> "$1"
}

log_error() {
	echo -e "[\033[1;31mERROR\033[m] \033[1;34m${FUNCNAME[1]}:\033[m $2" >&2
	echo "$(date -u '+%d-%m-%Y %T %Z') [ERROR] ${FUNCNAME[1]}: $2" >> "$1"
}
