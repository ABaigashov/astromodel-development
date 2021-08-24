if [[ ./renamer.sh != $BASH_SOURCE ]]; then
	echo "Oh no... You should run this script from scripts root: ./renamer.sh"
	exit 1
fi
cd ..
for dates in * ; do
	if [[ -d $dates && $dates != scripts ]]; then
		for object in $dates/* ; do
			counter=1001
			for photo in $object/* ; do
				if [[ fits == ${photo##*.} ]]; then
					mv -n $photo $object/${counter:1:3}.fits
					let counter++
				fi
				if [[ txt == ${photo##*.} ]]; then
					mv -n $photo $object/settings.txt
				fi
				if [[ csv == ${photo##*.} ]]; then
					mv -n $photo $object/histogram.csv
				fi
			done
		done
	fi
done