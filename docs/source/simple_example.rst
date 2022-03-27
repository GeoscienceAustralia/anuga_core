Simple Example
==============

Here we discuss the structure and operation of a
script called `runup.py` (which is available in the `examples`
directory of `anuga_core`.

This example carries out the solution of the shallow-water wave
equation in the simple case of a configuration comprising a flat
bed, sloping at a fixed angle in one direction and having a
constant depth across each line in the perpendicular direction.

The example demonstrates the basic ideas involved in setting up a
complex scenario. In general the user specifies the geometry
(bathymetry and topography), the initial water level, boundary
conditions such as tide, and any forcing terms that may drive the
system such as rainfall, abstraction of water, wind stress or atmospheric pressure gradients.
Frictional resistance from the different terrains in the model is
represented by predefined forcing terms. In this example, the
boundary is reflective on three sides and a time dependent wave on
one side.

The present example represents a simple scenario and does not
include any forcing terms, nor is the data taken from a file as it
would typically be.

The conserved quantities involved in the
problem are stage (absolute height of water surface),
$x$-momentum and $y$-momentum. Other quantities
involved in the computation are the friction and elevation.

Water depth can be obtained through the equation:

depth = stage - elevation

Outline of the Program
----------------------

In outline, `runup.py` performs the following steps:

   * Sets up a triangular mesh.
   * Sets certain parameters governing the mode of
         operation of the model, specifying, for instance,
         where to store the model output.
   * Inputs various quantities describing physical measurements, such
         as the elevation, to be specified at each mesh point (vertex).
   * Sets up the boundary conditions.
   * Carries out the evolution of the model through a series of time
         steps and outputs the results, providing a results file that can
         be viewed.


The Code
--------

For reference we include below the complete code listing for
`runup.py`. Subsequent paragraphs provide a
'commentary' that describes each step of the program and explains it
significance.


.. literalinclude:: ../../examples/simple_examples/runup.py

Establishing the Domain
-----------------------

The very first thing to do is import the various modules, of which the
`anuga` module is the most important:

.. code-block:: runup

    import anuga

Then we need to set up the triangular mesh to be used for the
scenario. This is carried out through the statement:

.. code-block:: runup

    domain = anuga.rectangular_cross_domain(10, 5, len1=10.0, len2=5.0)


The above assignment sets up a $10 \times
5$ rectangular mesh, triangulated in a regular way with boundary tags `left`, `right`,
         `top` or `bottom`.

It is also possible to set up a domain from ``first principles'' using \code{points}, \code{vertices} and \code{boundary} via the assignment:
\begin{verbatim}
domain = anuga.Domain(points, vertices, boundary)
\end{verbatim}
where
\begin{itemize}
   \item a list \code{points} giving the coordinates of each mesh point,
   \item a list \code{vertices} specifying the three vertices of each triangle, and
   \item a dictionary \code{boundary} that stores the edges on
         the boundary and associates with each a symbolic tag.
         The edges are represented as pairs (i, j) where i refers
         to the triangle id and j to the edge id of that triangle.
         Edge ids are enumerated from 0 to 2 based on the id of the vertex opposite.
\end{itemize}

(For more details on symbolic tags, see page
\pageref{ref:tagdescription}.)

An example of a general unstructured mesh and the associated data
structures \code{points}, \code{vertices} and \code{boundary} is
given in Section \ref{sec:meshexample}.



This creates an instance of the \class{Domain} class, which
represents the domain of the simulation. Specific options are set at
this point, including the basename for the output file and the
directory to be used for data:

\begin{verbatim}
domain.set_name('runup')
domain.set_datadir('.')
\end{verbatim}

In addition, the following statement could be used to state that
quantities \code{stage}, \code{xmomentum} and \code{ymomentum} are
to be stored at every timestep and \code{elevation} only once at
the beginning of the simulation:

\begin{verbatim}
domain.set_quantities_to_be_stored({
'stage': 2, 'xmomentum': 2, 'ymomentum': 2, 'elevation': 1})
\end{verbatim}

However, this is not necessary, as the above is the default behaviour.

\section{Initial Conditions}

The next task is to specify a number of quantities that we wish to
set for each mesh point. The class \class{Domain} has a method
\method{set\_quantity}, used to specify these quantities. It is a
flexible method that allows the user to set quantities in a variety
of ways -- using constants, functions, numeric arrays, expressions
involving other quantities, or arbitrary data points with associated
values, all of which can be passed as arguments. All quantities can
be initialised using \method{set\_quantity}. For a conserved
quantity (such as \code{stage, xmomentum, ymomentum}) this is called
an \emph{initial condition}. However, other quantities that aren't
updated by the equation are also assigned values using the same
interface. The code in the present example demonstrates a number of
forms in which we can invoke \method{set\_quantity}.

\subsection{Elevation}

The elevation, or height of the bed, is set using a function
defined through the statements below, which is specific to this
example and specifies a particularly simple initial configuration
for demonstration purposes:

\begin{verbatim}
def topography(x, y):
    return -x/2
\end{verbatim}

This simply associates an elevation with each point \code{(x, y)} of
the plane.  It specifies that the bed slopes linearly in the
\code{x} direction, with slope $-\frac{1}{2}$,  and is constant in
the \code{y} direction.

Once the function \function{topography} is specified, the quantity
\code{elevation} is assigned through the simple statement:

\begin{verbatim}
domain.set_quantity('elevation', topography)
\end{verbatim}

NOTE: If using function to set \code{elevation} it must be vector
compatible. For example, using square root will not work.

\subsection{Friction}

The assignment of the friction quantity (a forcing term)
demonstrates another way we can use \method{set\_quantity} to set
quantities -- namely, assign them to a constant numerical value:

\begin{verbatim}
domain.set_quantity('friction', 0.1)
\end{verbatim}

This specifies that the Manning friction coefficient is set to 0.1
at every mesh point.

\subsection{Stage}

The stage (the height of the water surface) is related to the
elevation and the depth at any time by the equation:

\begin{verbatim}
stage = elevation + depth
\end{verbatim}

For this example, we simply assign a constant value to \code{stage},
using the statement:

\begin{verbatim}
domain.set_quantity('stage', -0.4)
\end{verbatim}

which specifies that the surface level is set to a height of $-0.4$,
i.e.\ 0.4 units (metres) below the zero level.

Although it is not necessary for this example, it may be useful to
digress here and mention a variant to this requirement, which allows
us to illustrate another way to use \method{set\_quantity} -- namely,
incorporating an expression involving other quantities. Suppose,
instead of setting a constant value for the stage, we wished to
specify a constant value for the \emph{depth}. For such a case we
need to specify that \code{stage} is everywhere obtained by adding
that value to the value already specified for \code{elevation}. We
would do this by means of the statements:

\begin{verbatim}
h = 0.05    # Constant depth
domain.set_quantity('stage', expression='elevation + %f' % h)
\end{verbatim}

That is, the value of \code{stage} is set to $\code{h} = 0.05$ plus
the value of \code{elevation} already defined.

The reader will probably appreciate that this capability to
incorporate expressions into statements using \method{set\_quantity}
greatly expands its power. See Section \ref{sec:initial conditions} for more
details.

\section{Boundary Conditions}\index{boundary conditions}

The boundary conditions are specified as follows:

\begin{verbatim}
Br = anuga.Reflective_boundary(domain)
Bt = anuga.Transmissive_boundary(domain)
Bd = anuga.Dirichlet_boundary([0.2, 0.0, 0.0])
Bw = anuga.Time_boundary(domain=domain,
                   f=lambda t: [(0.1*sin(t*2*pi)-0.3)*exp(-t), 0.0, 0.0])
\end{verbatim}

The effect of these statements is to set up a selection of different
alternative boundary conditions and store them in variables that can be
assigned as needed. Each boundary condition specifies the
behaviour at a boundary in terms of the behaviour in neighbouring
elements. The boundary conditions introduced here may be briefly described as
follows:
\begin{itemize}
    \item \textbf{Reflective boundary}\label{def:reflective boundary}
          Returns same \code{stage} as in its neighbour volume but momentum
          vector reversed 180 degrees (reflected).
          Specific to the shallow water equation as it works with the
          momentum quantities assumed to be the second and third conserved
          quantities. A reflective boundary condition models a solid wall.
    \item \textbf{Transmissive boundary}\label{def:transmissive boundary}
          Returns same conserved quantities as
          those present in its neighbour volume. This is one way of modelling
          outflow from a domain, but it should be used with caution if flow is
          not steady state as replication of momentum at the boundary
          may cause numerical instabilities propagating into the domain and
          eventually causing \anuga to crash. If this occurs,
          consider using e.g.\ a Dirichlet boundary condition with a stage value
          less than the elevation at the boundary.
    \item \textbf{Dirichlet boundary}\label{def:dirichlet boundary} Specifies
          constant values for stage, $x$-momentum and $y$-momentum at the boundary.
    \item \textbf{Time boundary}\label{def:time boundary} Like a Dirichlet
          boundary but with behaviour varying with time.
\end{itemize}

\label{ref:tagdescription}Before describing how these boundary
conditions are assigned, we recall that a mesh is specified using
three variables \code{points}, \code{vertices} and \code{boundary}.
In the code we are discussing, these three variables are returned by
the function \code{rectangular}. The example given in
Section \ref{sec:realdataexample} illustrates another way of
assigning the values, by means of the function
\code{create_mesh_from_regions}.

These variables store the data determining the mesh as follows. (You
may find that the example given in Section \ref{sec:meshexample}
helps to clarify the following discussion, even though that example
is a \emph{non-rectangular} mesh.)
\begin{itemize}
    \item The variable \code{points} stores a list of 2-tuples giving the
          coordinates of the mesh points.
    \item The variable \code{vertices} stores a list of 3-tuples of
          numbers, representing vertices of triangles in the mesh. In this
          list, the triangle whose vertices are \code{points[i]},
          \code{points[j]}, \code{points[k]} is represented by the 3-tuple
          \code{(i, j, k)}.
    \item The variable \code{boundary} is a Python dictionary that
          not only stores the edges that make up the boundary but also assigns
          symbolic tags to these edges to distinguish different parts of the
          boundary. An edge with endpoints \code{points[i]} and
          \code{points[j]} is represented by the 2-tuple \code{(i, j)}. The
          keys for the dictionary are the 2-tuples \code{(i, j)} corresponding
          to boundary edges in the mesh, and the values are the tags are used
          to label them. In the present example, the value \code{boundary[(i, j)]}
          assigned to \code{(i, j)]} is one of the four tags
          \code{'left'}, \code{'right'}, \code{'top'} or \code{'bottom'},
          depending on whether the boundary edge represented by \code{(i, j)}
          occurs at the left, right, top or bottom of the rectangle bounding
          the mesh. The function \code{rectangular} automatically assigns
          these tags to the boundary edges when it generates the mesh.
\end{itemize}

The tags provide the means to assign different boundary conditions
to an edge depending on which part of the boundary it belongs to.
(In Section \ref{sec:realdataexample} we describe an example that
uses different boundary tags -- in general, the possible tags are entirely selectable by the user when generating the mesh and not
limited to 'left', 'right', 'top' and 'bottom' as in this example.)
All segments in bounding polygon must be tagged. If a tag is not supplied, the default tag name 'exterior' will be assigned by \anuga.

Using the boundary objects described above, we assign a boundary
condition to each part of the boundary by means of a statement like:

\begin{verbatim}
domain.set_boundary({'left': Br, 'right': Bw, 'top': Br, 'bottom': Br})
\end{verbatim}

It is critical that all tags are associated with a boundary condition in this statement.
If not the program will halt with a statement like:

\begin{verbatim}
Traceback (most recent call last):
  File "mesh_test.py", line 114, in ?
    domain.set_boundary({'west': Bi, 'east': Bo, 'north': Br, 'south': Br})
  File "X:\inundation\sandpits\onielsen\anuga_core\source\anuga\
    abstract_2d_finite_volumes\domain.py", line 505, in set_boundary
    raise msg
ERROR (domain.py): Tag "exterior" has not been bound to a boundary object.
All boundary tags defined in domain must appear in the supplied dictionary.
The tags are: ['ocean', 'east', 'north', 'exterior', 'south']
\end{verbatim}

The command \code{set_boundary} stipulates that, in the current example, the right
boundary varies with time, as defined by the lambda function, while the other
boundaries are all reflective.

The reader may wish to experiment by varying the choice of boundary
types for one or more of the boundaries. (In the case of \code{Bd}
and \code{Bw}, the three arguments in each case represent the
\code{stage}, $x$-momentum and $y$-momentum, respectively.)

\begin{verbatim}
Bw = anuga.Time_boundary(domain=domain,
f=lambda t: [(0.1*sin(t*2*pi)-0.3), 0.0, 0.0])
\end{verbatim}

\section{Evolution}\index{evolution}

The final statement:

\begin{verbatim}
for t in domain.evolve(yieldstep=0.1, duration=10.0):
    print domain.timestepping_statistics()
\end{verbatim}

causes the configuration of the domain to 'evolve', over a series of
steps indicated by the values of \code{yieldstep} and
\code{duration}, which can be altered as required.  The value of
\code{yieldstep} controls the time interval between successive model
outputs.  Behind the scenes more time steps are generally taken.

\section{Output}

The output is a NetCDF file with the extension \code{.sww}. It
contains stage and momentum information and can be used with the
\anuga viewer \code{anuga\_viewer} to generate a visual
display (see Section \ref{sec:anuga_viewer}). See Section \ref{sec:file formats}
(page \pageref{sec:file formats}) for more on NetCDF and other file
formats.

The following is a listing of the screen output seen by the user
when this example is run:

\verbatiminputunderscore{examples/runupoutput.txt}


\section{How to Run the Code}

The code can be run in various ways:
\begin{itemize}
  \item{from a Windows or Unix command line} as in\ \code{python runup.py}
  \item{within the Python IDLE environment}
  \item{within emacs}
  \item{within Windows, by double-clicking the \code{runup.py}
  file.}
\end{itemize}


\section{Exploring the Model Output}

The following figures are screenshots from the \anuga visualisation
tool \code{anuga_viewer}. Figure \ref{fig:runupstart} shows the domain
with water surface as specified by the initial condition, $t=0$.
Figure \ref{fig:runup2} shows later snapshots for $t=2.3$ and
$t=4$ where the system has been evolved and the wave is encroaching
on the previously dry bed.

\code{anuga_viewer} is described in more detail in Section \ref{sec:anuga_viewer}.

\begin{figure}[htp]
  \centerline{\includegraphics[width=75mm, height=75mm]
    {graphics/bedslopestart.jpg}}
  \caption{Runup example viewed with the \anuga viewer}
  \label{fig:runupstart}
\end{figure}

\begin{figure}[htp]
  \centerline{
    \includegraphics[width=75mm, height=75mm]{graphics/bedslopeduring.jpg}
    \includegraphics[width=75mm, height=75mm]{graphics/bedslopeend.jpg}
   }
  \caption{Runup example viewed with ANGUA viewer}
  \label{fig:runup2}
\end{figure}

\clearpage