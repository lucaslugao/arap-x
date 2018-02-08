from flask import Flask, send_from_directory, request
import os
import json
import numpy as np
from scipy import linalg
from scipy import sparse
from sksparse.cholmod import cholesky
from flask_sockets import Sockets
import datetime

print("ARAP")

app = Flask(__name__)
sockets = Sockets(app)

def analyseMatrix(M):
  print('Negative eigenvalues:', [(v,i) for (v,i) in enumerate(linalg.eigvals(M)) if v < 0])
  print('Symmetric: ', np.allclose(M, M.T))

def getAngles(vertices, triangle):
    """ Calculate the angles given a vertice list and a triangle of ids """
    a = vertices[triangle[0]]
    b = vertices[triangle[1]]
    c = vertices[triangle[2]]

    def angleBetween(a, b, c):
        ac = a - c
        bc = b - c
        return np.arccos(ac.dot(bc) / (linalg.norm(ac) * linalg.norm(bc)))

    return [angleBetween(b, c, a),
            angleBetween(c, a, b),
            angleBetween(a, b, c)]

def calculateWeigths(mesh):
  """ Calculate the weights ij based on the contagents """
  mesh['N'] = [set() for i in range(mesh['nV'])]
  mesh['weights'] = [dict() for i in range(mesh['nV'])]

  for f in mesh['faces']:
    mesh['N'][f[0]].add(f[1])
    mesh['N'][f[0]].add(f[2])

    mesh['N'][f[1]].add(f[0])
    mesh['N'][f[1]].add(f[2])

    mesh['N'][f[2]].add(f[0])
    mesh['N'][f[2]].add(f[1])

    angles = getAngles(mesh['p'], f)

    cot = [1/(2 * np.tan(angles[i])) for i in range(3)]

    mesh['weights'][f[0]][f[1]] = mesh['weights'][f[0]].get(f[1], 0) + cot[2]
    mesh['weights'][f[1]][f[0]] = mesh['weights'][f[1]].get(f[0], 0) + cot[2]

    mesh['weights'][f[0]][f[2]] = mesh['weights'][f[0]].get(f[2], 0) + cot[1]
    mesh['weights'][f[2]][f[0]] = mesh['weights'][f[2]].get(f[0], 0) + cot[1]

    mesh['weights'][f[1]][f[2]] = mesh['weights'][f[1]].get(f[2], 0) + cot[0]
    mesh['weights'][f[2]][f[1]] = mesh['weights'][f[2]].get(f[1], 0) + cot[0]

def calculateEnergy(mesh):
  """ Calculate the energy and energy derivative for debugging purposes """
  E = 0
  for i in range(mesh['nV']): 
    for j in mesh['N'][i]:
      E += mesh['weights'][i][j] * linalg.norm((mesh['p_prime'][i] - mesh['p_prime'][j]) - mesh['R'][i].dot(mesh['p'][i] - mesh['p'][j]))**2
  D = np.array([0.0,0.0,0.0])
  for i in range(mesh['nV']):
    if i not in mesh['fixed']: 
      for j in mesh['N'][i]:
        D += 4 * mesh['weights'][i][j] * ((mesh['p_prime'][i] - mesh['p_prime'][j]) - (mesh['R'][i]+mesh['R'][j]).dot(mesh['p'][i] - mesh['p'][j])/2)
  return E, linalg.norm(D)

def calculateRigidRotations(mesh):
  """ Calculate the rigid rotations that minimises the mesh energy """
  mesh['R'] = [None for i in range(mesh['nV'])]

  # sign flipping matrix
  correction = np.eye(3)
  correction[2][2] = -1

  for i in range(mesh['nV']):
    P = np.zeros((3, len( mesh['N'][i])))
    P_prime = np.zeros((3, len( mesh['N'][i])))
    D = np.zeros((len( mesh['N'][i]), len( mesh['N'][i])))
    for c, j in enumerate(mesh['N'][i]):
      D[c][c] = mesh['weights'][i][j]
      for r in range(3):
        P[r][c] = mesh['p'][j][r]
        P_prime[r][c] = mesh['p_prime'][j][r]

    S = np.matmul(np.matmul(P, D), P_prime.T)
    U, s, V = linalg.svd(S)
    V = V.T
    mesh['R'][i] = np.matmul(V, U.T)
    if linalg.det(mesh['R'][i]) < 0:
        mesh['R'][i] = np.matmul(np.matmul(V, correction), U.T)
    
          
def prefactorLMatrix(mesh):
  """ Prefactor the L sparse matrix """
  row = []
  col = []
  data = []
  for i in range(mesh['nV']):
      row.append(i)
      col.append(i)
      if i in mesh['fixed']:
          data.append(1)
      else:
          data.append(sum([mesh['weights'][i][j] for j in mesh['N'][i]]))
          for j in mesh['N'][i]:
            # Don't include a fixed point in the LHS of the system
            if j not in mesh['fixed']:
              row.append(i)
              col.append(j)
              data.append(-mesh['weights'][i][j])

  L = sparse.csc_matrix((data, (row, col)), shape=(mesh['nV'], mesh['nV']))
  mesh['factor'] = cholesky(L)
  mesh['L'] = L
  
def deform(mesh):
  """ Deform the vertices positions using the prefactore L matrix """
  def get_b_entry(i):
    """ Auxiliary function to calculate the i-th row of b """
    if i in mesh['fixed']:
      return mesh['p_prime'][i]
    bi = np.zeros(3)
    for j in mesh['N'][i]:
      wij = mesh['weights'][i][j]
      Rij = mesh['R'][i] + mesh['R'][j]
      pij = mesh['p'][i] - mesh['p'][j]
      bi += (wij/2) * Rij.dot(pij)
      # Adds the correction of the fixed point to b
      if j in mesh['fixed']:
        bi += mesh['weights'][i][j] * mesh['p_prime'][j]
    return bi
  b = np.array([ get_b_entry(i) for i in range(mesh['nV'])])
  mesh['p_prime'] = mesh['factor'](b)


@sockets.route('/engine')
def engine(ws):
    """ WebSocket endpoint """
    print("Client connected {}".format(ws))
    mesh = {}
    while not ws.closed:
      msg = ws.receive()
      if msg is None:
        break
      msgObj = json.loads(msg)
      now = datetime.datetime.utcnow().isoformat() + 'Z'
      print(now, msgObj['action'])
      
      # Reads and process each action accordingly 

      if(msgObj['action'] == 'create'):
        mesh['p'] = [np.array(x, dtype= np.dtype('Float64')) for x in msgObj['vertices']] 
        mesh['faces'] = msgObj['faces']
        mesh['nV'] = len(mesh['p'])

      if(msgObj['action'] == 'update' or msgObj['action'] == 'create'):
        mesh['fixed'] = set(msgObj['fixed'])
        calculateWeigths(mesh)
        prefactorLMatrix(mesh)

      if(msgObj['action'] == 'deform'):
        mesh['p_prime'] = [np.array(x, dtype= np.dtype('Float64')) for x in msgObj['vertices']]
        iterations = msgObj['iterations'] if 'iterations' in msgObj else 4
        energy = [0,0]
        def printEnergy(label, energy):
          newEnergy =  calculateEnergy(mesh)
          head = 'Energy [{}] = '.format(label)
          print(head, newEnergy,'\n' + ''.join([' ' for i in range(len(head))]) ,'Delta = ', newEnergy[0] - energy[0])
          return newEnergy
        for i in range(iterations): 
          calculateRigidRotations(mesh)
          energy = printEnergy('rot {}'.format(i), energy)
          deform(mesh)
          energy = printEnergy('def {}'.format(i), energy)
          new_p_prime = [{'x': float(v[0]), 'y': float(v[1]), 'z': float(v[2])} for v in mesh['p_prime']]
          ws.send(json.dumps(new_p_prime))

    print("Client disconnected {}".format(ws))


@app.route('/static/<path:path>')
def serve_static(path):
  """ Serve the static files in the /static/ http endpoint """
  return send_from_directory('static', path)

@app.route('/', defaults={'file': 'index.html'})
def serve_results(file):
    """ Special route for the index.html file """
    return send_from_directory('static', file)

# Server configuration and binding
from gevent import pywsgi
from geventwebsocket.handler import WebSocketHandler

host = "0.0.0.0"
port = int(os.environ['PORT']) if 'PORT' in os.environ.keys() else 80
print("Trying to server at {}:{}".format(host, port))
server = pywsgi.WSGIServer((host, port), app, handler_class=WebSocketHandler)
server.serve_forever()