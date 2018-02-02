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

def getAngles(vertices, triangle):
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

  if False:
    for i in range(mesh['nV']):
      for j in mesh['N'][i]:
        assert(mesh['weights'][i][j] == mesh['weights'][j][i])
        assert(i in mesh['N'][j])

def calculateEnergy(mesh):
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
  mesh['R'] = [np.zeros((3, 3)) for i in range(mesh['nV'])]

  correction = np.eye(3)
  correction[2][2] = -1
  for i in range(mesh['nV']):
      S = np.zeros((3, 3))
      for j in mesh['N'][i]:
          S += mesh['weights'][i][j] * (mesh['p'][i] - mesh['p'][j]) * \
              np.transpose([mesh['p_prime'][i] - mesh['p_prime'][j]])
    
      U, s, V = linalg.svd(S)
      mesh['R'][i] = np.matmul(V.T, U.T)
      if linalg.det(mesh['R'][i]) < 0:
          mesh['R'][i] = np.matmul(np.matmul(V.T, correction), U.T)

def prefactorLMatrix(mesh):
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
            if j not in mesh['fixed']:
              row.append(i)
              col.append(j)
              data.append(-mesh['weights'][i][j])

  L = sparse.csc_matrix((data, (row, col)), shape=(mesh['nV'], mesh['nV']))
  mesh['factor'] = cholesky(L)
  mesh['L'] = L
  

def deform(mesh):
  def get_b_entry(i):
    if i in mesh['fixed']:
      return mesh['p_prime'][i]
    bi = np.zeros(3)
    for j in mesh['N'][i]:
      wij = mesh['weights'][i][j]
      Rij = mesh['R'][i] + mesh['R'][j]
      pij = mesh['p'][i] - mesh['p'][j]
      bi += (wij/2) * Rij.dot(pij)
      if j in mesh['fixed']:
        bi += mesh['weights'][i][j] * mesh['p_prime'][j]
    return bi
  b = np.array([ get_b_entry(i) for i in range(mesh['nV'])])

          
  #new_p = linalg.solve(mesh['L'].todense(), b)
  new_p = mesh['factor'](b)
  
  #for j in mesh['fixed']:
  #  print(linalg.norm(new_p[j]- mesh['p_prime'][j]))
  if False:
    #Lp = mesh['L'].todense().dot(new_p)
    for i in range(mesh['nV']):
      if i not in mesh['fixed']:
        lhs = np.zeros(3)
        rhs = np.zeros(3)
        for j in mesh['N'][i]:
          lhs += mesh['weights'][i][j] * new_p[i]
          lhs -= mesh['weights'][i][j] * new_p[j]
          rhs += (mesh['weights'][i][j]/2) * (mesh['R'][i] + mesh['R'][j]).dot(mesh['p'][i] - mesh['p'][j])
        #print(linalg.norm(lhs-Lp[i]) < 1.0e-9, lhs,'<->', Lp[i], ' => ' )
        #print(lhs,'<->', rhs, ' => ', linalg.norm(lhs-rhs))
  mesh['p_prime'] = new_p


@sockets.route('/engine')
def engine(ws):
    print("Client connected {}".format(ws))
    mesh = {}
    while not ws.closed:
      msg = ws.receive()
      if msg is None:
        break
      msgObj = json.loads(msg)
      now = datetime.datetime.utcnow().isoformat() + 'Z'
      print(now, msgObj['action'])
      
      if(msgObj['action'] == 'create'):

        mesh['p'] = [np.array(x) for x in msgObj['vertices']] 
        mesh['faces'] = msgObj['faces']
        mesh['fixed'] = set(msgObj['fixed'])
        mesh['nV'] = len(mesh['p'])

        calculateWeigths(mesh)
        prefactorLMatrix(mesh)
      else:
        mesh['p_prime'] = [np.array(x) for x in msgObj['vertices']]
        iterations = msgObj['iterations'] if 'iterations' in msgObj else 4
        for i in range(iterations): 
          calculateRigidRotations(mesh)
          if i == 0:
            print('Energy [r] = ', calculateEnergy(mesh))
          deform(mesh)
          print('Energy [', i, '] = ', calculateEnergy(mesh))
          ws.send(json.dumps([{'x': float(v[0]), 'y': float(v[1]), 'z': float(v[2])} for v in mesh['p_prime']]))

    print("Client disconnected {}".format(ws))

@app.route('/static/<path:path>')
def serve_static(path):
  return send_from_directory('static', path)


@app.route('/', defaults={'file': 'index.html'})
def serve_results(file):
    return send_from_directory('static', file)

from gevent import pywsgi
from geventwebsocket.handler import WebSocketHandler

host = "0.0.0.0"
port = int(os.environ['PORT']) if 'PORT' in os.environ.keys() else 80
print("Trying to server at {}:{}".format(host, port))
server = pywsgi.WSGIServer((host, port), app, handler_class=WebSocketHandler)
server.serve_forever()