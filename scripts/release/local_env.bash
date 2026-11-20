#!/bin/bash
# Node-local wrappers for the build's Python and R environments.
# Fontconfig, matplotlib and XDG caches live under LOCALDIR so renders never
# lock caches on shared network storage. Versions are checked against
# environment/*-lock.txt by scripts/release/check_environment.py.
# Usage: local_env.bash LOCALDIR PYTHON_ENV_PREFIX R_ENV_PREFIX
#   then make PYTHON=LOCALDIR/bin/python RSCRIPT=LOCALDIR/bin/Rscript ...
set -euo pipefail
L=$1; PY_ENV=$2; R_ENV=$3
mkdir -p "$L/bin" "$L/xdg_cache" "$L/mpl"
for pair in "py:$PY_ENV" "r:$R_ENV"; do
  tag=${pair%%:*}; env=${pair#*:}
  d="$L/fc_$tag"; mkdir -p "$d/cache"
  cp -p "$env/var/cache/fontconfig/"*.cache-9 "$d/cache/" 2>/dev/null || true
  sed -e "s#<cachedir>$env/var/cache/fontconfig</cachedir>#<cachedir>$d/cache</cachedir>#" \
      -e "s#<cachedir>~/.fontconfig</cachedir>##" \
      -e "s#<include ignore_missing=\"yes\">conf.d</include>#<include ignore_missing=\"yes\">$env/etc/fonts/conf.d</include>#" \
      "$env/etc/fonts/fonts.conf" > "$d/fonts.conf"
  grep -q "$d/cache" "$d/fonts.conf"
  ! grep -q "$env/var/cache" "$d/fonts.conf"
done
cat > "$L/bin/python" <<W
#!/bin/bash
export FONTCONFIG_FILE=$L/fc_py/fonts.conf XDG_CACHE_HOME=$L/xdg_cache MPLCONFIGDIR=$L/mpl
exec $PY_ENV/bin/python "\$@"
W
cat > "$L/bin/Rscript" <<W
#!/bin/bash
export FONTCONFIG_FILE=$L/fc_r/fonts.conf XDG_CACHE_HOME=$L/xdg_cache
export LD_LIBRARY_PATH=$R_ENV/lib:\${LD_LIBRARY_PATH:-}
exec $R_ENV/bin/Rscript "\$@"
W
chmod +x "$L/bin/python" "$L/bin/Rscript"
