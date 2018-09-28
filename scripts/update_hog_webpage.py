
from os import listdir
from os.path import isfile, join


path_webiste='./website/'

# https://stackoverflow.com/questions/3207219/how-do-i-list-all-files-of-a-directory
onlyfiles = [f for f in listdir(mypath) if isfile(join(path_webiste, f)) and 'phobius' not in f ]


header='<a href="../index.html">Home</a><br><br>'

search_box='\n\n\t\t<div>\n\t\t\t<form name="f" action="http://lmse.utoronto.ca/aybrah/hl.php" method="post">\n\t\t\t\t<b>Lookup String:</b> <input type="text" value="" size="15" name=text>\n\t\t\t\t<input type="submit" value="Search" name=search>\n\t\t\t</form>\n\t\t</div>\n\n'


for file_html in onlyfiles:
	page=open(path_webiste+file_html,'r').read()
	page=page.replace('#banner { background: #bdbdbd; height:200px;','#banner { background: #bdbdbd; height:225px;')
	page=page.replace(header,header+search_box)
	open(path_webiste+file_html,'w').write(page)