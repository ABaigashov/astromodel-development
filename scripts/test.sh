result=${PWD##*/}
selfname=${0##*/}

if [[ "astro-model" != "$result" ]]; then
  echo "Oh no... You should run this script from repository root: ./scripts/$selfname [ARGS]"
  exit 1
fi

echo $result
