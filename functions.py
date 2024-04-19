def H(A, q):
  return as_tensor([[0, 1], [1190476.19048 * A ** 0.5 - (q / A) ** 2, 2 * q / A]])

def F(A, q):
  return as_vector([q, 793650.79365 * A ** 1.5 - 6.34921 + q ** 2 / A])

def B(A, q):
  return as_vector([0, k_r * q / A])

def S(A, q):
  return as_vector([0, k_r * q / A])

def dS_dU(A, q):
  return as_tensor([[0, 0], [- k_r * q / A ** 2, k_r / A]])

def U(A, q):
  return np.array([A, q])

def H_vec(A, q):
  return np.array([[0, 1], [1190476.19048 * A ** 0.5 - (q / A) ** 2, 2 * q / A]])

def B_vec(A, q):
  return np.array([0, k_r * q / A])

def S_vec(A, q):
  return np.array([0, k_r * q / A])

def c_alpha(A, q, alpha = 1):
  return (1190476.19048 * A ** 0.5 + (q / A) ** 2 * alpha * (alpha - 1)) ** 0.5

def l1(A, q, alpha = 1):
  return np.array([c_alpha(A, q, alpha) - alpha * q / A, 1.])

def l2(A, q, alpha = 1):
  return np.array([- c_alpha(A, q, alpha) - alpha * q / A, 1.])

def CC(A, q, u_der, x):
  du_dz = project(u_der, V_der)(x)
  return U(A, q) - dt * H_vec(A, q) @ du_dz - dt * B_vec(A, q) + dt * f([x, x, x])

def inlet_bc(A, q, u_der):
  q_inlet = (np.dot(l2(A(x_left), q(x_left)), CC(A(x_left), q(x_left), u_der, x_left)) - l2(A(x_left), q(x_left))[0] * U(A(x_left), q(x_left))[0] ) / l2(A(x_left), q(x_left))[1]
  
  return q_inlet

def outlet_bc(A, q, u_der):
  matrix = np.array([[1, -1], [c_alpha(A(x_right), q(x_right), alpha) + alpha * q(x_right) / A(x_right), c_alpha(A(x_right), q(x_right), alpha) - alpha * q(x_right) / A(x_right)]])
  array = np.array([np.dot(l1(A(x_right), q(x_right)), CC(A(x_right), q(x_right), u_der, x_right)), np.dot(l2(A(x_right), q(x_right)), U(A(x_right), q(x_right)) - dt * S_vec(A(x_right), q(x_right)))])

  A_out, q_out = matrix @ array / (2 * c_alpha(A(x_right), q(x_right)))

  return A_out, q_out

def update_time_step(t, dt):
  t += dt
  q_exact.t = t
  A_exact.t = t
  dU_dt.t = t
  S.t = t
  dF_dx.t = t
  dU_dt_dt.t = t
  dS_dt.t = t
  dF_dx_dt.t = t

  return t