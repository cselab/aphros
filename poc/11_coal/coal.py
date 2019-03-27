def xmf(path, hfn, field, ncx, ncy, ncz, ncc):
	extent = 1.
	 
	nx = ncx + 1
	ny = ncy + 1
	nz = ncz + 1
	bx = 0. 
	by = 0. 
	bz = 0. 
	dx = extent / max(ncx, ncy, ncz)
	dy = dx
	dz = dx
	prec = 8
	hfs = "data"
	
	t = '''<?xml version='1.0' ?>
<!DOCTYPE Xdmf SYSTEM 'Xdmf.dtd' []>
<Xdmf Version='2.0'>
 <Domain>
   <Grid GridType='Uniform'>
     <Topology TopologyType='3DCORECTMesh' Dimensions='{nz} {ny} {nx}'/>
     <Geometry GeometryType='ORIGIN_DXDYDZ'>
       <DataItem Name='Origin' Dimensions='3' NumberType='Float' Precision='8' Format='XML'>
        {bz} {by} {bx}
       </DataItem>
       <DataItem Name='Spacing' Dimensions='3' NumberType='Float' Precision='8' Format='XML'>
        {dz} {dy} {dx}
       </DataItem>
     </Geometry>
     <Attribute Name='{field}' AttributeType='Scalar' Center='Cell'>
       <DataItem Dimensions='{ncz} {ncy} {ncx} {ncc}' NumberType='Float' Precision='{prec}' Format='HDF'>
        ./{hfn}:/{hfs}
       </DataItem>
     </Attribute>
   </Grid>
 </Domain>
</Xdmf>'''.format(**locals())
	open(path, 'w').write(t)
