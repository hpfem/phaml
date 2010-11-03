def plot_mesh_mpl(polygons=None, polynomial_orders=None, edges_only=False):
    """
    Plots an hp-FEM mesh (optionally including the element orders).

    The ``polygons`` is a dictionary mapping element ids to an array of (x, y)
    coordinates of the polygonial boundary of the element (3 points for a
    triangle, 4 points for a quad).

    The ``polynomial_orders`` is a dictionary mapping element ids to their
    polynomial degrees.

    The function returns a Matplotlib ``Figure`` instance, that you can then
    use to save it to a file, or do some other things with it.

    Example:

    >>> from numpy import array
    >>> from femhub.plot import plot_mesh_mpl
    >>> f = plot_mesh_mpl({
            0: array([[0, 0], [1, 0], [1, 1]]),
            1: array([[1, 1], [0.5, 0.5], [0.5, 1]]),
            }, {0: 1, 1: 4})
    >>> f.savefig("a.png")

    """
    from matplotlib import pyplot
    from matplotlib.path import Path
    from matplotlib.patches import PathPatch
    from matplotlib.patches import Rectangle
    colors_old = {
            1: '#000684',
            2: '#3250fc',
            3: '#36c4ee',
            4: '#04eabc',
            5: '#62ff2a',
            6: '#fdff07',
            7: '#ffa044',
            8: '#ff1111',
            9: '#b02c2c',
            10: '#820f97',
            }
    colors = {
            0: '#7f7f7f',
            1: '#7f2aff',
            2: '#2a2aff',
            3: '#2a7fff',
            4: '#00d4aa',
            5: '#00aa44',
            6: '#abc837',
            7: '#ffd42a',
            8: '#c87137',
            9: '#c83737',
            10: '#ff0000',
            }
    fig = pyplot.figure()
    sp = fig.add_subplot(111)
    nodes = []
    for el_id in polygons:
        x = list(polygons[el_id][:, 0])
        y = list(polygons[el_id][:, 1])
        nodes.extend(zip(x, y))
        x.append(x[0])
        y.append(y[0])
        vertices = zip(x, y)
        codes = [Path.MOVETO] + [Path.LINETO]*(len(vertices)-2) + \
                    [Path.CLOSEPOLY]
        p = Path(vertices, codes)
        if edges_only:
            color = "white"
            linewidth = 2
        else:
            if polynomial_orders is None:
                color = colors[0]
            else:
                color = colors.get(polynomial_orders[el_id], "#ffff00")
            linewidth = 1
        patch = PathPatch(p, facecolor=color, lw=linewidth,
                edgecolor='#000000')
        sp.add_patch(patch)
    show_legend = polynomial_orders is not None

    if show_legend:
        # Create legend
        def split_nodes():
            x = []
            y = []

            if isinstance(nodes, dict):
                _nodes = nodes.items()
            else:
                _nodes = enumerate(nodes)
            for k, pnt in _nodes:
                x.append(pnt[0])
                y.append(pnt[1])

            return (x, y)

        def get_max(what='x'):
            x, y = split_nodes()

            if what == 'x':
                return max(x)
            else:
                return max(y)

        def get_min(what='x'):
            x, y = split_nodes()

            if what == 'x':
                return min(x)
            else:
                return min(y)

        maxX = get_max('x')
        maxY = get_max('y')

        minX = get_min('x')
        minY = get_min('y')

        dy = (maxY - minY) / 20
        dx = (maxX - minX) / 20

        y = minY + dy
        x = maxX + dx

        ord = polynomial_orders.items()
        order_list = []
        for k,v in ord:
            order_list.append(v)
        m = max(order_list)

        for k,c in colors.items():
            if k <= m :
                p = Rectangle(xy=(x,y), width=dx, height=dy, fill=True, facecolor=c)
                sp.add_patch(p)
                sp.text(x + dx + (dx/2), y + (dy/4), str(k))
                y += dy
            else:
                break

        sp.text(x, y + (dy/2), str('Orders'))
    sp.set_title("Mesh")
    sp.set_aspect("equal")
    sp.autoscale_view()
    return sp.figure

def plot_sln_mayavi(x, y, mesh, sln_values, colorbar=False):
    """
    Plot a solution using mayavi.

    Example:

    >>> from numpy import array
    >>> from femhub.plot import plot_sln_mayavi
    >>> f = plot_sln_mayavi([0, 1, 1], [0, 0, 1], [1, 2, 3])
    >>> f.savefig("a.png")

    """
    from enthought.mayavi import mlab
    #mlab.options.offscreen = True
    mlab.clf()
    #mlab.options.show_scalar_bar = False
    z = [0] * len(x)
    mlab.triangular_mesh(x, y, z, mesh, scalars=sln_values)
    engine = mlab.get_engine()
    image = engine.current_scene
    image.scene.background = (1.0, 1.0, 1.0)
    image.scene.foreground = (0.0, 0.0, 0.0)
    if colorbar:
        mlab.colorbar(orientation="vertical")
    mlab.view(0, 0)
    return mlab

def convert_mesh(x, y, elems, elems_orders):
    """
    Convert the mesh from Phaml representation to femhub representation.
    """
    from numpy import array
    polygons = {}
    for n, elem in enumerate(elems):
        polygons[n] = array([ [x[i-1], y[i-1]] for i in elem ])
    orders = {}
    for n, order in enumerate(elems_orders):
        orders[n] = order
    return polygons, orders
