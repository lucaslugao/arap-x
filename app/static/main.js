var container, stats;
var camera, scene, renderer;
var handlerObjects = [];
var positions = [];
var options = {iterations: 4};

var fixed_handle_toggle = true;
var hovering = null;

var handlerGeometry = new THREE.CubeGeometry(3, 3, 3);
var handlerMaterials = [
  new THREE.MeshLambertMaterial({color: 0xcc0000}),
  new THREE.MeshLambertMaterial({color: 0x00ff00})
];
var transformControl;

var mainObject = null;

if (window.location.protocol == "https:") {
  var ws_scheme = "wss://";
} else {
  var ws_scheme = "ws://"
};

var ws = new WebSocket(ws_scheme + location.host + "/engine");

ws.addEventListener('open', function() {
  console.log('Python backend connected!');
  init();
  animate();
});



function loadOBJ(path, cb) {
  var loader = new THREE.OBJLoader();
  loader.load(
      path,
      function(object) {
        object.traverse(function(child) {
          if (child instanceof THREE.Mesh) {
            cb(child);
          }
        });

      },
      function(xhr) {
        console.log((xhr.loaded / xhr.total * 100) + '% loaded');
      },
      function(error) {
        console.log('An error happened');
      });
}

function addHandler(target, vertexID) {
  var material = handlerMaterials[0];
  var object = new THREE.Mesh(handlerGeometry, material);

  object.position.x = target.geometry.vertices[vertexID].x;
  object.position.y = target.geometry.vertices[vertexID].y;
  object.position.z = target.geometry.vertices[vertexID].z;
  object.state = 0;
  object.dirty = true;

  object.toggleState = function() {
    object.state = 1 - object.state;
    object.material = handlerMaterials[object.state];
    if(mainObject !== null)
      mainObject.dirty = true;
  };

  if (object.position.y < 10) object.toggleState();

  object.vertexID = vertexID;
  object.move = function() {
    target.geometry.vertices[vertexID].set(
        object.position.x, object.position.y, object.position.z);
    target.geometry.verticesNeedUpdate = true;
  };

  object.fix = false;
  object.castShadow = true;
  object.receiveShadow = true;
  scene.add(object);
  handlerObjects.push(object);
}

function packageVertices(vertices) {
  return vertices.map(function(v) {
    return [v.x, v.y, v.z];
  });
}

function packageFaces(faces) {
  return faces.map(function(v) {
    return [v.a, v.b, v.c];
  });
}

function getFixed(handlers) {
  var fixed = handlers
  .filter(function(v) {
    return v.state == 1;
  })
  .map(function(v) {
    return v.vertexID;
  });
  
  if(fixed.length == 0)
    fixed.push(0);
  
  return fixed;
}
/*

handlerObjects[183].position.set(10,235,255); handlerObjects[183].move();
handlerObjects[183].toggleState(); requestDeformation();

*/
function load(scene) {
  var url_string = window.location.href
  var url = new URL(url_string);
  var c = url.searchParams.get("model");
  if(c === null)
    c = 'cube'
  loadOBJ('static/models/' + c  +'.obj', function(object) {
    object.material = new THREE.MeshPhongMaterial({
      color: 0xffffff,
      specular: 0xffffff,
      shininess: 20,
      morphTargets: true,
      vertexColors: THREE.FaceColors,
      flatShading: true
    });
    /*object.material =
        new THREE.MeshLambertMaterial({color: 0x0000ff, wireframe: true});*/
    object.material.side = THREE.DoubleSide;
    object.geometry = new THREE.Geometry().fromBufferGeometry(object.geometry);
    object.geometry.mergeVertices();
    scene.add(object);
    object.geometry.vertices
        .filter(function(v, i) {
          return true;
        })
        .forEach(function(v, i) {
          addHandler(object, i);
        });
    mainObject = object;
    object.geometry.originalPackedVertices =
        packageVertices(object.geometry.vertices)

    ws.addEventListener('message', function(m) {
      console.log('[ws:onmessage] Received new vertices!');
      vertices = JSON.parse(m.data);
      vertices.forEach(function(v, i) {
        mainObject.geometry.vertices[i].set(v.x, v.y, v.z);
        handlerObjects[i].position.set(v.x, v.y, v.z);
      });
      mainObject.geometry.verticesNeedUpdate = true;
    })
    sendObject();
  });
}

function sendObject() {
  if (mainObject !== null) {
    mainObject.dirty = false;
    console.log('[sendObject] New object sent!');
    ws.send(JSON.stringify({
      action: 'create',
      vertices: mainObject.geometry.originalPackedVertices,
      faces: packageFaces(mainObject.geometry.faces),
      fixed: getFixed(handlerObjects)
    }));
  }
}


function updateObject() {
  if (mainObject !== null) {
    mainObject.dirty = false;
    console.log('[updateObject] New object sent!');
    ws.send(JSON.stringify({
      action: 'update',
      fixed: getFixed(handlerObjects)
    }));
  }
}


function requestDeformation(iterations = options.iterations) {
  if (mainObject !== null) {
    if (mainObject.dirty) updateObject();
    console.log('[requestDeformation] New fixed positions sent!');
    ws.send(JSON.stringify({
      action: 'deform',
      vertices: packageVertices(mainObject.geometry.vertices),
      iterations: iterations
    }));
  }
}

function init() {
  container = document.getElementById('container');

  scene = new THREE.Scene();
  scene.background = new THREE.Color(0xf0f0f0);

  camera = new THREE.PerspectiveCamera(
      70, window.innerWidth / window.innerHeight, 1, 10000);
  camera.position.set(0, 250, 1000);
  scene.add(camera);

  hemiLight = new THREE.HemisphereLight(0xffffff, 0xffffff, 0.6);
  hemiLight.color.setHSL(0.6, 1, 0.6);
  hemiLight.groundColor.setHSL(0.095, 1, 0.75);
  hemiLight.position.set(0, 50, 0);
  scene.add(hemiLight);
  hemiLightHelper = new THREE.HemisphereLightHelper(hemiLight, 10);
  scene.add(hemiLightHelper);

  dirLight = new THREE.DirectionalLight(0xffffff, 1);
  dirLight.color.setHSL(0.1, 1, 0.95);
  dirLight.position.set(-1, 1.75, 1);
  dirLight.position.multiplyScalar(30);
  scene.add(dirLight);
  dirLight.castShadow = true;
  dirLight.shadow.mapSize.width = 2048;
  dirLight.shadow.mapSize.height = 2048;
  var d = 50;
  dirLight.shadow.camera.left = -d;
  dirLight.shadow.camera.right = d;
  dirLight.shadow.camera.top = d;
  dirLight.shadow.camera.bottom = -d;
  dirLight.shadow.camera.far = 3500;
  dirLight.shadow.bias = -0.0001;
  dirLightHeper = new THREE.DirectionalLightHelper(dirLight, 10);
  scene.add(dirLightHeper);

  var planeGeometry = new THREE.PlaneGeometry(2000, 2000);
  planeGeometry.rotateX(-Math.PI / 2);
  var planeMaterial = new THREE.ShadowMaterial({opacity: 0.2});

  var plane = new THREE.Mesh(planeGeometry, planeMaterial);
  plane.position.y = 0;
  plane.receiveShadow = true;
  scene.add(plane);

  var helper = new THREE.GridHelper(2000, 100);
  helper.position.y = 1;
  helper.material.opacity = 0.25;
  helper.material.transparent = true;
  scene.add(helper);

  renderer = new THREE.WebGLRenderer({antialias: true});
  renderer.setPixelRatio(window.devicePixelRatio);
  renderer.setSize(window.innerWidth, window.innerHeight);
  renderer.shadowMap.enabled = true;
  container.appendChild(renderer.domElement);

  stats = new Stats();
  container.appendChild(stats.dom);

  var gui = new dat.GUI();

  gui.add({Calculate: requestDeformation}, 'Calculate');
  gui.add({Toggle_non_fixed: function(){
    fixed_handle_toggle = !fixed_handle_toggle;
    handlerObjects.forEach(function(v){
      if(v.state == 0)
        v.visible = fixed_handle_toggle;
    })
  }}, 'Toggle_non_fixed');
  gui.add(options, 'iterations', 1, 35).step(1);
  gui.open();


  var controls = new THREE.OrbitControls(camera, renderer.domElement);
  controls.damping = 0.2;
  controls.addEventListener('change', render);

  controls.addEventListener('start', function() {
    cancelHideTransorm();
  });

  controls.addEventListener('end', function() {
    delayHideTransform();
  });

  transformControl = new THREE.TransformControls(camera, renderer.domElement);
  transformControl.addEventListener('change', render);
  scene.add(transformControl);

  // Hiding transform situation is a little in a mess :()
  transformControl.addEventListener('change', function(e) {
    cancelHideTransorm();
  });

  transformControl.addEventListener('mouseDown', function(e) {
    cancelHideTransorm();
  });

  transformControl.addEventListener('mouseUp', function(e) {
    delayHideTransform();
    /*if(e.target.object.state == 1)
      requestDeformation();*/
    
  });

  transformControl.addEventListener('objectChange', function(e) {
    e.target.object.move();
  });

  var dragcontrols =
      new THREE.DragControls(handlerObjects, camera, renderer.domElement);
  dragcontrols.enabled = false;

  dragcontrols.addEventListener('hoveron', function(event) {
    transformControl.attach(event.object);
    cancelHideTransorm();
    hovering = event.object;
  });

  dragcontrols.addEventListener('hoveroff', function(event) {
    delayHideTransform();
    hovering = null;
  });


  document.body.onkeyup = function(e) {
    if (e.keyCode == 32) {
      requestDeformation();
    } else if (e.keyCode == 65) {
      if (hovering !== null) {
        hovering.toggleState();
        /*if(hovering.state == 1)
          requestDeformation();*/
      }
    } else if (e.keyCode == 82) {
      if (hovering !== null) {
        var p = mainObject.geometry.originalPackedVertices[hovering.vertexID];
        mainObject.geometry.vertices[hovering.vertexID].set(p[0], p[1], p[2]);
        hovering.position.set(p[0], p[1], p[2]);
        mainObject.geometry.verticesNeedUpdate = true;
      }
    } else if (e.keyCode == 73) {
      if (hovering !== null) {
        console.log({
          id: hovering.vertexID,
          originalPos:
              mainObject.geometry.originalPackedVertices[hovering.vertexID],
          pos: mainObject.geometry.vertices[hovering.vertexID]
        });
      }
    }
  };

  var hiding;

  function delayHideTransform() {
    cancelHideTransorm();
    hideTransform();
  }

  function hideTransform() {
    hiding = setTimeout(function() {
      transformControl.detach(transformControl.object);
    }, 2500);
  }

  function cancelHideTransorm() {
    if (hiding) clearTimeout(hiding);
  }
  load(scene);
}

function animate() {
  requestAnimationFrame(animate);
  render();
  stats.update();
  transformControl.update();
}

function render() {
  renderer.render(scene, camera);
}