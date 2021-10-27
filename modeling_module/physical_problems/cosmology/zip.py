import os
import shutil
import zipfile
import matplotlib.pyplot as plt

x = [1, 3, 5]
y = [4, 5, 1]
z = [9, 1, 9]

try:
    os.mkdir('tmp')
except:
    print('The tmp folder is already exist')

for i in range(5):
    plt.plot(x, y)
    plt.savefig('tmp/pic1')
    plt.plot(y, x)
    plt.savefig('tmp/pic2')
    plt.plot(z, x)
    plt.savefig('tmp/pic3')

fantasy_zip = zipfile.ZipFile('archive.zip', 'w')

for folder, subfolders, files in os.walk('tmp'):

    for file in files:
        if file.endswith('.png'):
            fantasy_zip.write(os.path.join(folder, file), os.path.relpath(os.path.join(folder,file), 'tmp'), compress_type = zipfile.ZIP_DEFLATED)

fantasy_zip.close()
shutil.rmtree('tmp')
