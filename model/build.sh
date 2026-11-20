#!/usr/bin/env bash
# Build ICM-Hydro with gfortran into model/bin/icm_hydro.
set -e
cd "$(dirname "$0")"
OBJ=build/obj; mkdir -p "$OBJ" bin
FLAGS="-ffixed-line-length-none -ffree-line-length-none -O2 -fallow-argument-mismatch -std=legacy -w -J $OBJ"
todo=$(cat srclist.txt)
# a file that uses a module compiles only after that module: repeat passes until all compile
for pass in 1 2 3 4 5 6 7 8; do
  left=""
  for f in $todo; do
    o="$OBJ/$(echo "${f%.*}" | tr / _).o"
    gfortran -c $FLAGS "src/$f" -o "$o" 2> build/last_error.txt || left="$left $f"
  done
  n=$(echo $left | wc -w | tr -d ' ')
  echo "pass $pass: $n files left"
  [ "$n" -eq 0 ] && break
  [ "$(echo $left)" = "$(echo $todo)" ] && { echo "build stopped, see model/build/last_error.txt:$left"; exit 1; }
  todo=$left
done
gfortran -o bin/icm_hydro "$OBJ"/*.o
echo "built model/bin/icm_hydro"
