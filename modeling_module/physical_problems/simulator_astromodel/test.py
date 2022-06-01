
import json

f = open('init_file.json', )
data = json.load(f)
data2 = data["obj_1"]
print(data2["mass"])