function getShader(gl, id) {
  var shaderScript = document.getElementById(id);
  if (!shaderScript) {
    return null;
  }

  var str = "";
  var k = shaderScript.firstChild;
  while (k) {
    if (k.nodeType == 3) {
      str += k.textContent;
    }
    k = k.nextSibling;
  }

  var shader;
  if (shaderScript.type == "x-shader/x-fragment") {
    shader = gl.createShader(gl.FRAGMENT_SHADER);
  } else if (shaderScript.type == "x-shader/x-vertex") {
    shader = gl.createShader(gl.VERTEX_SHADER);
  } else {
    return null;
  }

  gl.shaderSource(shader, str);
  gl.compileShader(shader);

  if (!gl.getShaderParameter(shader, gl.COMPILE_STATUS)) {
    alert(gl.getShaderInfoLog(shader));
    return null;
  }

  return shader;
}

var mvMatrix;
function loadIdentity() { mvMatrix = Matrix.I(4); }
function multMatrix(m) { mvMatrix = mvMatrix.x(m); }
function mvTranslate(v) {
  var m = Matrix.Translation($V([v[0], v[1], v[2]])).ensure4x4();
  multMatrix(m);
}
var pMatrix;
function perspective(fovy, aspect, znear, zfar) {
  pMatrix = makePerspective(fovy, aspect, znear, zfar);
}

function ortho(left, right, bottom, top, znear, zfar){ 
	pMatrix = makeOrtho(left, right, bottom, top, znear, zfar)
}


function setMatrixUniforms() {
  gl.uniformMatrix4fv(shaderProgram.pMatrixUniform, false, new Float32Array(pMatrix.flatten()));
  gl.uniformMatrix4fv(shaderProgram.mvMatrixUniform, false, new Float32Array(mvMatrix.flatten()));
}




//implementation specific
// 
// 
// var triangleVertexPositionBuffer;
//   var triangleVertexColorBuffer;
//   var squareVertexPositionBuffer;
//   var squareVertexColorBuffer;
//   function initBuffers() {
//     triangleVertexPositionBuffer = gl.createBuffer();
//     gl.bindBuffer(gl.ARRAY_BUFFER, triangleVertexPositionBuffer);
//     var vertices = [
//          0.0,  1.0,  0.0,
//         -1.0, -1.0,  0.0,
//          1.0, -1.0,  0.0
//     ];
//     gl.bufferData(gl.ARRAY_BUFFER, new Float32Array(vertices), gl.STATIC_DRAW);
//     triangleVertexPositionBuffer.itemSize = 3;
//     triangleVertexPositionBuffer.numItems = 3;
// 
//     triangleVertexColorBuffer = gl.createBuffer();
//     gl.bindBuffer(gl.ARRAY_BUFFER, triangleVertexColorBuffer);
//     var colors = [
//         1.0, 0.0, 0.0, 1.0,
//         0.0, 1.0, 0.0, 1.0,
//         0.0, 0.0, 1.0, 1.0,
//     ];
//     gl.bufferData(gl.ARRAY_BUFFER, new Float32Array(colors), gl.STATIC_DRAW);
//     triangleVertexColorBuffer.itemSize = 4;
//     triangleVertexColorBuffer.numItems = 3;
// 
// 
//     squareVertexPositionBuffer = gl.createBuffer();
//     gl.bindBuffer(gl.ARRAY_BUFFER, squareVertexPositionBuffer);
//     vertices = [
//          1.0,  1.0,  0.0,
//         -1.0,  1.0,  0.0,
//          1.0, -1.0,  0.0,
//         -1.0, -1.0,  0.0
//     ];
//     gl.bufferData(gl.ARRAY_BUFFER, new Float32Array(vertices), gl.STATIC_DRAW);
//     squareVertexPositionBuffer.itemSize = 3;
//     squareVertexPositionBuffer.numItems = 4;
// 
//     squareVertexColorBuffer = gl.createBuffer();
//     gl.bindBuffer(gl.ARRAY_BUFFER, squareVertexColorBuffer);
//     colors = []
//     for (var i=0; i < 4; i++) {
//       colors = colors.concat([0.5, 0.5, 1.0, 1.0]);
//     }
//     gl.bufferData(gl.ARRAY_BUFFER, new Float32Array(colors), gl.STATIC_DRAW);
//     squareVertexColorBuffer.itemSize = 4;
//     squareVertexColorBuffer.numItems = 4;
//   }
// 
// 
// 
//   function drawScene() {
//     gl.viewport(0, 0, gl.viewportWidth, gl.viewportHeight);
//     gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);
// 
//     perspective(45, gl.viewportWidth / gl.viewportHeight, 0.1, 100.0);
//     loadIdentity();
// 
//     mvTranslate([-1.5, 0.0, -7.0])
// 
//     gl.bindBuffer(gl.ARRAY_BUFFER, triangleVertexPositionBuffer);
//     gl.vertexAttribPointer(shaderProgram.vertexPositionAttribute, triangleVertexPositionBuffer.itemSize, gl.FLOAT, false, 0, 0);
// 
//     gl.bindBuffer(gl.ARRAY_BUFFER, triangleVertexColorBuffer);
//     gl.vertexAttribPointer(shaderProgram.vertexColorAttribute, triangleVertexColorBuffer.itemSize, gl.FLOAT, false, 0, 0);
// 
//     setMatrixUniforms();
//     gl.drawArrays(gl.TRIANGLES, 0, triangleVertexPositionBuffer.numItems);
// 
//     mvTranslate([3.0, 0.0, 0.0])
//     gl.bindBuffer(gl.ARRAY_BUFFER, squareVertexPositionBuffer);
//     gl.vertexAttribPointer(shaderProgram.vertexPositionAttribute, squareVertexPositionBuffer.itemSize, gl.FLOAT, false, 0, 0);
// 
//     gl.bindBuffer(gl.ARRAY_BUFFER, squareVertexColorBuffer);
//     gl.vertexAttribPointer(shaderProgram.vertexColorAttribute, squareVertexColorBuffer.itemSize, gl.FLOAT, false, 0, 0);
// 
//     setMatrixUniforms();
//     gl.drawArrays(gl.TRIANGLE_STRIP, 0, squareVertexPositionBuffer.numItems);
//   }
