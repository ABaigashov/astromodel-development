from fenics import *
import numpy as np
from mshr import *
import numpy as np
import matplotlib.pyplot as plt
import sympy as sym
import os


# from mpl_toolkits.mplot3d import *
# from matplotlib import cm
# import matplotlib.animation as animation
############################
#####УРАВНЕНИЕ ПУАССОНА#####
############################
path = 'modeling_module/physical_problems/poisson_problem/'

class Task_maker():

    def __init__(self, config):
        self.config = config

        self.name_vector = []
        self.name_scalar = []
        self.description_vector = []
        self.description_scalar = []
        self.type_vector = []
        self.type_scalar = []
        self.vector_condition = []
        self.scalar_condition = []
        self.notation_scalar = []
        self.source = "0"
        self.kappa = "1"

        if self.config.mesh_name:
             self.mesh = path + 'mesh/' + self.config.mesh_name

        if self.config.source:
            self.source = self.config.source

        # if self.config.kappa:
        #     self.kappa = self.config.kappa

        # if self.config.vector_functions:
        #     self.conditions = self.config.vector_functions
        #     for boundaries in self.conditions:
        #         if boundaries.boundary_name:
        #             self.name_vector.append(boundaries.boundary_name)
        #         if boundaries.description:
        #             self.description_vector.append(boundaries.description)
        #         if boundaries.condition_type:
        #             self.type_vector.append(boundaries.condition_type)
        #         if boundaries.condition:
        #             self.vector_condition.append(boundaries.condition)

        if self.config.scalar_functions:
            self.conditions = self.config.scalar_functions
            for boundaries in self.conditions:
                if boundaries.name:
                    self.name_scalar.append(boundaries.name)
                if boundaries.description:
                    self.description_scalar.append(boundaries.description)
                if boundaries.condition_type:
                    self.type_scalar.append(boundaries.condition_type)
                if boundaries.condition:
                    self.scalar_condition.append(boundaries.condition)
                if boundaries.notation:
                    self.notation_scalar.append(boundaries.notation)
                else:
                    self.notation_scalar.append("EXPR")


class BVP_solver():

    def __init__(self, task, output):

        self.mesh = Mesh(task.mesh)
        self.output = output

        self.V = FunctionSpace(self.mesh, 'P', 1)
        #объявление искомых функций и пробных функций. Они являются частью V.
        self.u = Function(self.V)
        self.v = TestFunction(self.V)

        self.Dc = []
        self.Nc = []
        self.bx = []

        self.boundary_parts = MeshFunction('size_t', self.mesh, self.mesh.topology().dim() - 1)

        ds = Measure('ds', domain = self.mesh, subdomain_data = self.boundary_parts)

        k=0
        m=0
        for i in task.name_scalar:
            if task.type_scalar[k]=="Dirichlet":
                if task.notation_scalar[k]=="SYM":
                    u_D = sym.sympify(task.scalar_condition[k])
                    x, y, z, t = sym.symbols('x[0], x[1], x[2], t')
                    sub = [('x',x),('y',y),('z',z),('t',t)]
                    u_D = u_D.subs(sub)
                    u_code = sym.printing.ccode(u_D)
                    u_D = Expression(u_code, degree=2)
                else:
                    u_D = Expression(task.scalar_condition[k], degree=2)
                DC = DirichletBC(self.V, u_D, task.description_scalar[k])
                self.Dc.append(DC)
            if task.type_scalar[k]=="Neumann":
                if task.notation_scalar[k]=="SYM":
                    g_N = sym.sympify(task.scalar_condition[k])
                    x, y, z, t = sym.symbols('x[0], x[1], x[2], t')
                    sub = [('x',x),('y',y),('z',z),('t',t)]
                    g_N = g_N.subs(sub)
                    g_code = sym.printing.ccode(g_N)
                    g = Expression(g_code, degree=2)
                else:
                    g = Expression(task.scalar_condition[k], degree=2)
                bx1 = CompiledSubDomain(task.description_scalar[k])
                self.bx.append(bx1)
                self.bx[m].mark(self.boundary_parts, m)
                self.Nc.append(g*self.v*ds(m))
                m = m+1
            k = k+1

            #источник в правой части уравнения Пуассона
        self.f = Expression(task.source, degree=2, u=self.u)

        self.kappa1 = Expression(task.kappa, degree=2,u=self.u)


    def Solving_eq(self):
        #постановка вариационной задачи и ее решение с граничными условиями
        F = self.kappa1*dot(grad(self.u), grad(self.v))*dx - self.f*self.v*dx - sum(self.Nc)
        solve(F == 0, self.u, self.Dc)

        # Save solution to file in VTK format
        vtkfile = File(path+'results/'+'solution.pvd')
        vtkfile << self.u

        return f'{self.output}/solution.pvd'

class PointsPotential:

    def __init__(self, config, output, job):
        self.config = config
        self.output = output
        os.mkdir(self.output)
        self.job = job

    def start(self):

        self.parameters = self.config.mass_object
        #-------------Create DOLPHIN mesh and define function space----------
        ae = 1.496 * 10**11
        G = 6.67408 * 10**(-11)
        tol = 1E-14
        self.edge = 100
        self.ncells = 256

        domain =  Rectangle(Point(-self.edge, -self.edge),
                            Point(self.edge, self.edge))

        # generate mesh and determination of the Function Space
        mesh = generate_mesh(domain, self.ncells)

        P1 = FiniteElement('CG', triangle, 1)

        V = FunctionSpace(mesh, P1)
        Q = FunctionSpace(mesh, 'DG', 0)

        ex_L = ''
        ex_R = ''
        ex_H = ''
        ex_B = ''

        for p in self.parameters:
        	if p.y >= 0 and p.x >= 0:
        		ex_L += f'- G * {p.mass_obj} / pow(pow(x[1]-{p.y}, 2) + pow(-{self.edge}-{p.x}, 2), 0.5)'
        		ex_R += f'- G * {p.mass_obj} / pow(pow(x[1]-{p.y}, 2) + pow({self.edge}-{p.x}, 2), 0.5)'
        		ex_H += f'- G * {p.mass_obj} / pow(pow(x[0]-{p.x}, 2) + pow({self.edge}-{p.y}, 2), 0.5)'
        		ex_B += f'- G * {p.mass_obj} / pow(pow(x[0]-{p.x}, 2) + pow(-{self.edge}-{p.y}, 2), 0.5)'

        	elif y_obj[i] < 0 and x_obj[i] >= 0:
        		ex_L += f'- G * {p.mass_obj} / pow(pow(x[1]+{-p.y}, 2) + pow(-{self.edge}-{p.x}, 2), 0.5)'
        		ex_R += f'- G * {p.mass_obj} / pow(pow(x[1]+{-p.y}, 2) + pow({self.edge}-{p.x}, 2), 0.5)'
        		ex_H += f'- G * {p.mass_obj} / pow(pow(x[0]-{p.x}, 2) + pow({self.edge}+{-p.y}, 2), 0.5)'
        		ex_B += f'- G * {p.mass_obj} / pow(pow(x[0]-{p.x}, 2) + pow(-{self.edge}+{-p.y}, 2), 0.5)'

        	elif y_obj[i] >= 0 and x_obj[i] < 0:
        		ex_L += f'- G * {p.mass_obj} / pow(pow(x[1]-{p.y}, 2) + pow(-{self.edge}+{-p.x}, 2), 0.5)'
        		ex_R += f'- G * {p.mass_obj} / pow(pow(x[1]-{p.y}, 2) + pow({self.edge}+{-p.x}, 2), 0.5)'
        		ex_H += f'- G * {p.mass_obj} / pow(pow(x[0]+{-p.x}, 2) + pow({self.edge}-{p.y}, 2), 0.5)'
        		ex_B += f'- G * {p.mass_obj} / pow(pow(x[0]+{-p.x}, 2) + pow(-{self.edge}-{p.y}, 2), 0.5)'

        	elif y_obj[i] < 0 and x_obj[i] < 0:
        		ex_L += f'- G * {p.mass_obj} / pow(pow(x[1]+{-p.y}, 2) + pow(-{self.edge}+{-p.x}, 2), 0.5)'
        		ex_R += f'- G * {p.mass_obj} / pow(pow(x[1]+{-p.y}, 2) + pow({self.edge}+{-xp.x}, 2), 0.5)'
        		ex_H += f'- G * {p.mass_obj} / pow(pow(x[0]+{-p.x}, 2) + pow({self.edge}+{-p.y}, 2), 0.5)'
        		ex_B += f'- G * {p.mass_obj} / pow(pow(x[0]+{-p.x}, 2) + pow(-{self.edge}+{-p.y}, 2), 0.5)'

        # -------------Define boundary condition------------
        phi_L = Expression(ex_L, G=G, degree=2)
        phi_R = Expression(ex_R, G=G, degree=2)
        phi_H = Expression(ex_H, G=G, degree=2)
        phi_B = Expression(ex_B, G=G, degree=2)

        def boundary_L(x, on_boundary):
        	return on_boundary and near(x[0], -self.edge, tol)
        def boundary_R(x, on_boundary):
        	return on_boundary and near(x[0], self.edge, tol)
        def boundary_H(x, on_boundary):
        	return on_boundary and near(x[1], self.edge, tol)
        def boundary_B(x, on_boundary):
        	return on_boundary and near(x[1], -self.edge, tol)

        bc_L = DirichletBC(V, phi_L, boundary_L)
        bc_R = DirichletBC(V, phi_R, boundary_R)
        bc_H = DirichletBC(V, phi_H, boundary_H)
        bc_B = DirichletBC(V, phi_B, boundary_B)

        bc = [bc_L, bc_R, bc_H, bc_B]

        # -------------Define variational problem------------
        phi = Function(V)
        v = TestFunction(V)

        rho = []
        rho_sum = Expression('0', degree=2)
        for p in self.parameters:
        	condition = f'pow(pow(x[0] - {p.x}, 2) + pow(x[1] - {p.y}, 2), 0.5) <= {p.r} ? 3 * {p.mass_obj} / (4 * pi * pow({p.r}, 3)): 0'
        	rho_i = Expression(condition, degree=2, pi=np.pi)
        	rho.append(rho_i)
        	rho_sum += rho_i


        J = Expression('pow(x[0]*x[0] + x[1]*x[1], 0.5)', degree=2)

        Func = (J*phi.dx(0)*v.dx(0) + J*phi.dx(1)*v.dx(1) - J*4*np.pi*rho_sum*v)*dx

        #---------------------Compute solution---------------------
        solve(Func == 0, phi, bc)

        #--------------------------Ploting------------------------
        #get array componets and triangulation :
        v = phi.compute_vertex_values(mesh)
        x = mesh.coordinates()[:,0]
        y = mesh.coordinates()[:,1]
        t = mesh.cells()

        ax = plt.axes()

        cm = plt.get_cmap('viridis')
        c = ax.tricontourf(x, y, t, v, 10, cmap=cm)
        p = ax.triplot(x, y, t, '-', color='k', lw=0.2, alpha=0.4)

        plt.xlim(-self.edge, self.edge)
        plt.ylim(-self.edge, self.edge)
        plt.axis('equal')

        #-----------Output in the file-------------------
        plt.savefig(f"{self.output}/potential_func.png")

        return f"{self.output}/potential_func.png"



#пример задания граничных условий через привычные переменные x,y на языке Python
# x, y = sym.symbols('x[0], x[1]')
# u_D = 1 + x + 2*y
#
# #пример сложного выражения для правой части уравнения Пуассона в терминах sympy
# f = - sym.diff(q(u_D)*sym.diff(u_D, x), x) - sym.diff(q(u_D)*sym.diff(u_D, y), y)
# f = sym.simplify(f)
#
# #конвертация в код C++
# u_code = sym.printing.ccode(u_D)
# f_code = sym.printing.ccode(f)
# print('u =', u_code)
# print('f =', f_code)
#
# #граничные условия Дирихле
# u_D = Expression(u_code, degree=2)
# def boundary(x, on_boundary):
#     return on_boundary
# bc = DirichletBC(V, u_D, boundary)
#
# f = Expression(f_code, degree=1)
# F = q(u)*dot(grad(u), grad(v))*dx - f*v*dx
#
# solve(F == 0, u, bc)
#
# # Plot solution and mesh
# plot(u)
# plot(mesh)
#
# # Save solution to file in VTK format
# vtkfile = File('ex5/solution.pvd')
# vtkfile << u
