f = open('files.out','r')
fout = open('convert_AAS.exe','w')

for line in f:
  line = line.rsplit()
  args = line[0].split('.')
  lineout = 'convert '+line[0]+' '+args[0]+'.png \n'
  fout.write(lineout)

fout.close()
f.close()


convert -delay 10 -loop 0 ls -v *.png animation.gif
