import xml.dom.minidom as dom

d = dom.Document()
e = d.createElement('svg')
e.setAttribute('xmlns', 'http://www.w3.org/2000/svg')
e.setAttribute('version', '1.1')
e.setAttribute('viewBox', '0 0 100 100')
e.setAttribute('width', '100%')
e.setAttribute('height', '50px')
e.setAttribute('style', 'width:100%; height:100%;')
e.setAttribute('preserveAspectRatio', 'none')

r = d.createElement('rect')
e.appendChild(r)
a = {'x':0, 'y':0, 'width':100, 'height':100, 'fill':'red'}
for k,v in a.items(): r.setAttribute(k,str(v))

#d.appendChild(e)
print e.toprettyxml()
