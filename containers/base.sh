#!/usr/bin/env sh

# Adds paths and useful functions to environment
# Should work for dash and bash

if [ -d ${ENV_DIR} ]; then
  for i in ${ENV_DIR}/*.sh; do
    if [ -r $i ]; then
      . $i
    fi
  done
  unset i
fi

# Adds an environment variable to env.
add_env () {
  echo "export ${1}=\"${2}\"" >> "${ENV_FILE}"
}

# Adds a path to the environment, without overwriting it.
prepend_path () {
  echo "export ${1}=\"${2}:\${${1}}\"" >> "${ENV_FILE}"
  eval export ${1}="${2}:\${${1}}"
}

# Adds some packages to be installed by apt
add_runtime_dep () {
  for dep in ${@}
  do
    echo ${dep} >> "${APT_REQUIREMENTS_FILE}"
  done
}

apt_install_from_file () {
  for f in ${@}
  do
    xargs -a "${f}" -r -- apt-get install -y --no-install-recommends
  done
}
