import base64, json
ui=open('src/ui.html',encoding='utf8').read(); eng=open('src/engine.js',encoding='utf8').read()
b64=lambda p: base64.b64encode(open(p,'rb').read()).decode()
S={
 'cad':{'kind':'raster','name':'cad_plan.png','src':'github.com/grebtsew/FloorplanToBlender3d (Images/Examples/example.png)','data':'data:image/png;base64,'+b64('samples/cad_plan.png'),
        'preset':{'pW':13.5,'pN':2,'pH':3.0,'pL':5.0,'pO':3.5,'pR':'0'}},
 'sketch':{'kind':'raster','name':'sketch.jpg','src':'github.com/grebtsew/FloorplanToBlender3d (Images/Examples/example6.png)','data':'data:image/jpeg;base64,'+b64('samples/sketch.jpg'),
        'preset':{'pW':6.0,'pN':1,'pH':3.0,'pL':5.0,'pO':3.5,'pR':'0'}},
 'dxf':{'kind':'dxf','name':'Test.dxf','src':'github.com/omarsamy3/DXF_To_Revit (CADToRevit/DXF Files/Test.dxf)','data':open('samples/test.dxf',encoding='utf8',errors='replace').read(),
        'preset':{'pU':'auto','pT':0.25,'pN':1,'pH':6.0,'pL':8.0,'pO':6.0,'pR':'0'}}}
html=ui.replace('/*__ENGINE__*/',eng).replace('/*__SAMPLES__*/',json.dumps(S,ensure_ascii=False).replace('</','<\\/'))
open('index.html','w',encoding='utf8').write(html); print(len(html))
