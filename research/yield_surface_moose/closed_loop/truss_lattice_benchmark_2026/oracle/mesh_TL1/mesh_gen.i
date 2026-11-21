[Mesh]
  allow_renumbering = true
  [e0]
    type = ElementGenerator
    nodal_positions = '0.0 0.0 0  1.0 0.0 0'
    element_connectivity = '0 1'
    elem_type = EDGE2
    subdomain_id = 1
  []
  [e1]
    type = ElementGenerator
    input = e0
    nodal_positions = '0.0 0.0 0  0.5 0.8660254037844386 0'
    element_connectivity = '0 1'
    elem_type = EDGE2
    subdomain_id = 2
  []
  [e2]
    type = ElementGenerator
    input = e1
    nodal_positions = '0.0 0.0 0  -0.5 0.8660254037844386 0'
    element_connectivity = '0 1'
    elem_type = EDGE2
    subdomain_id = 3
  []
  [e3]
    type = ElementGenerator
    input = e2
    nodal_positions = '1.0 0.0 0  2.0 0.0 0'
    element_connectivity = '0 1'
    elem_type = EDGE2
    subdomain_id = 1
  []
  [e4]
    type = ElementGenerator
    input = e3
    nodal_positions = '1.0 0.0 0  1.5 0.8660254037844386 0'
    element_connectivity = '0 1'
    elem_type = EDGE2
    subdomain_id = 5
  []
  [e5]
    type = ElementGenerator
    input = e4
    nodal_positions = '1.0 0.0 0  0.5 0.8660254037844386 0'
    element_connectivity = '0 1'
    elem_type = EDGE2
    subdomain_id = 6
  []
  [e6]
    type = ElementGenerator
    input = e5
    nodal_positions = '0.5 0.8660254037844386 0  1.5 0.8660254037844386 0'
    element_connectivity = '0 1'
    elem_type = EDGE2
    subdomain_id = 4
  []
  [e7]
    type = ElementGenerator
    input = e6
    nodal_positions = '0.5 0.8660254037844386 0  1.0 1.7320508075688772 0'
    element_connectivity = '0 1'
    elem_type = EDGE2
    subdomain_id = 2
  []
  [e8]
    type = ElementGenerator
    input = e7
    nodal_positions = '0.5 0.8660254037844386 0  0.0 1.7320508075688772 0'
    element_connectivity = '0 1'
    elem_type = EDGE2
    subdomain_id = 6
  []
  [e9]
    type = ElementGenerator
    input = e8
    nodal_positions = '1.5 0.8660254037844386 0  2.5 0.8660254037844386 0'
    element_connectivity = '0 1'
    elem_type = EDGE2
    subdomain_id = 4
  []
  [e10]
    type = ElementGenerator
    input = e9
    nodal_positions = '1.5 0.8660254037844386 0  2.0 1.7320508075688772 0'
    element_connectivity = '0 1'
    elem_type = EDGE2
    subdomain_id = 5
  []
  [e11]
    type = ElementGenerator
    input = e10
    nodal_positions = '1.5 0.8660254037844386 0  1.0 1.7320508075688772 0'
    element_connectivity = '0 1'
    elem_type = EDGE2
    subdomain_id = 3
  []
  [e12]
    type = ElementGenerator
    input = e11
    nodal_positions = '-1.0 -1.0 0  -2.0 -1.0 0'
    element_connectivity = '0 1'
    elem_type = EDGE2
    subdomain_id = 99
  []
  [merge]
    type = MeshRepairGenerator
    input = e12
    fix_node_overlap = true
  []
  [ns_pin]
    type = ExtraNodesetGenerator
    input = merge
    new_boundary = 'pin'
    coord = '0.0 0.0 0'
  []
  [ns_c1]
    type = ExtraNodesetGenerator
    input = ns_pin
    new_boundary = 'c1'
    coord = '-1.0 -1.0 0'
  []
  [ns_c2]
    type = ExtraNodesetGenerator
    input = ns_c1
    new_boundary = 'c2'
    coord = '-2.0 -1.0 0'
  []
[]
