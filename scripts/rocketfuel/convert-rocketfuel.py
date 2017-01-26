import argparse
import re

import networkx as nx
import matplotlib.pyplot as plt
import xml.etree.cElementTree as et

from xml.dom import minidom

def convert_to_scn(topology, req_size=15):

    root = et.Element("topology")

    x = 0
    y = 0
    for node in topology.nodes():

        # for each node create a sub-element of type 'router'
        node_id = topology.node[node]['node_id']
        domain = et.SubElement(root, "router", id = node_id, tier = node_id.split(".")[0], index = node_id.split(".")[1])

        # within each node, add links
        for neigh_node in topology.neighbors(node):

            edge_id = str(topology[node][neigh_node]['edge_nr'])
            link = et.SubElement(domain, "link", local = edge_id, remote = edge_id, rrouter = topology.node[neigh_node]['node_id'])

            for i in xrange(req_size):
                fwd_size_dist = et.SubElement(link, "fwd_size_dist", size = str(i + 1)).text = "0.00"

        x = x + 1
        y = y + 1

    xmlstr = minidom.parseString(et.tostring(root)).toprettyxml(indent="    ")
    with open("topology.scn", "w") as f:
        f.write(xmlstr)

def parse_pop_level_isp_map(rocketfuel_path):

    topology = nx.Graph()
    comment_char = '#'
    # extract the isp prefix (i.e. '<isp-number>:')
    isp_nr = rocketfuel_path.split("/")[-1] + ":"

    # extract the node names and edges
    node_id = 0
    for line in open(rocketfuel_path.rstrip("/") + "/edges", "r").readlines():

        if comment_char in line:
            # split on comment char, keep only the part before
            line, _ = line.split(comment_char, 1)
            line = line.strip()

        if len(line) == 0:
            continue

        try:

            # extract the edge number
            edge_nr = int(re.findall(" \d+", line)[-1])
            print("edge_nr = %d" % (edge_nr))

            # extract head and tail nodes of the edge
            pop_head = line.split(" -> ")[0]
            pop_tail = line.split(" -> ")[1]
            pop_tail = re.findall("[1-9]\d*:\D*, \S*", pop_tail)[0]
            print("pop_head = %s, pop_tail = %s" % (pop_head, pop_tail))

            # add nodes to pop-level topology graph

            if pop_head not in topology:
                topology.add_node(pop_head, node_id = str(node_id) + ".0")
                node_id = node_id + 1

            if pop_tail not in topology:
                topology.add_node(pop_tail, node_id = str(node_id) + ".0")
                node_id = node_id + 1

            # add edge
            topology.add_edge(pop_head, pop_tail, edge_nr = edge_nr, edge_latency = 0.0, edge_weight = 0.0)
            print("topology[%s][%s] = %d" % (pop_head, pop_tail, topology[pop_head][pop_tail]['edge_nr']))

        except IndexError:
            raise ValueError('Invalid input file. Parsing failed '\
                             'while trying to parse an internal node')

    # extract the latency from edges.lat
    for line in open(rocketfuel_path.rstrip("/") + "/edges.lat", "r").readlines():

        if comment_char in line:
            # split on comment char, keep only the part before
            line, _ = line.split(comment_char, 1)
            line = line.strip()

        if len(line) == 0:
            continue

        try:

            # latency is a decimal number at the end of the line
            print("line = %s" % (line))
            print("match = %s" % (re.findall("\d+\.*\d*", line)[-1]))
            edge_latency = float(re.findall("\d+\.*\d*", line)[-1])
            print("edge_latency = %f" % (edge_latency))
            # extract head and tail nodes of the edge
            pop_head = line.split(" -> ")[0]
            pop_tail = line.split(" -> ")[1]
            pop_tail = re.findall("[1-9]\d*:\D*, \S*", pop_tail)[0]

            # update the edge attributes
            topology[pop_head][pop_tail]['edge_latency'] = edge_latency
            print("topology[%s][%s] = %f" % (pop_head, pop_tail, topology[pop_head][pop_tail]['edge_latency']))

        except IndexError:
            raise ValueError('Invalid input file. Parsing failed '\
                             'while trying to parse an internal node')

    # extract the weight (whatever that is...) from edges.wt
    for line in open(rocketfuel_path.rstrip("/") + "/edges.wt", "r").readlines():

        if comment_char in line:
            # split on comment char, keep only the part before
            line, _ = line.split(comment_char, 1)
            line = line.strip()

        if len(line) == 0:
            continue

        try:

            # weight is also a decimal number, at the end of the line
            edge_weight = float(re.findall("\d+\.*\d*", line)[-1])
            print("edge_weight = %f" % (edge_weight))
            # extract head and tail nodes of the edge
            pop_head = line.split(" -> ")[0]
            pop_tail = line.split(" -> ")[1]
            pop_tail = re.findall("[1-9]\d*:\D*, \S*", pop_tail)[0]

            # update the edge attributes
            topology[pop_head][pop_tail]['edge_weight'] = edge_weight
            print("topology[%s][%s] = %f" % (pop_head, pop_tail, topology[pop_head][pop_tail]['edge_weight']))

        except IndexError:
            raise ValueError('Invalid input file. Parsing failed '\
                             'while trying to parse an internal node')

    return topology

# parser for RocketFuel ISP (router-level) maps .cch files
# adapted from http://fnss.github.io/ (email: fnss.dev@gmail.com) 
def parse_isp_map(rocketfuel_path):
    """
    Parse a network topology from RocketFuel ISP map file.

    The ASes provided by the RocketFuel dataset are the following:

    +------+---------------------+-------+--------+------------+------------+
    | ASN  | Name                | Span  | Region | Nodes (r1) | Nodes (r0) |
    +======+=====================+=======+========+============+============+
    | 1221 | Telstra (Australia) | world | AUS    |  2999      |  378 (318) |
    | 1239 | Sprintlink (US)     | world | US     |  8352      |  700 (604) |
    | 1755 | EBONE (Europe)      | world | Europe |   609      |  172       |
    | 2914 | Verio (US)          | world | US     |  7109      | 1013       |
    | 3257 | Tiscali (Europe)    | world | Europe |   855      |  248 (240) |
    | 3356 | Level 3 (US)        | world | US     |  3447      |  652       |
    | 3967 | Exodus (US)         | world | US     |   917      |  215 (201) |
    | 4755 | VSNL (India)        | world | India  |   121      |   12       |
    | 6461 | Abovenet (US)       | world | US     |     0      |  202       |
    | 7018 | AT&T (US)           | world | US     | 10152      |  656 (631) |
    +------+---------------------+-------+--------+------------+------------+

    Parameters
    ----------
    rocketfuel_path : str
        The path of the file containing the RocketFuel map. It should have
        extension .cch

    Returns
    -------
    topology : DirectedTopology
        The object containing the parsed topology.

    Notes
    -----
    The returned topology is always directed. If an undirected topology is
    desired, convert it using the DirectedTopology.to_undirected() method.

    Each node of the returned graph has the following attributes:
     * **type**: string
     * **location**: string (optional)
     * **address**: string
     * **r**: int
     * **backbone**: boolean (optional)

    Each edge of the returned graph has the following attributes:
     * type : string, which can either be *internal* or *external*

    If the topology contains self-loops (links starting and ending in the same
    node) they are stripped from the topology.

    Raises
    ------
    ValueError
        If the provided file cannot be parsed correctly.

    Examples
    --------
    >>> import fnss
    >>> topology = fnss.parse_rocketfuel_isp_map('1221.r0.cch')
    """

    topology = nx.Graph()
    comment_char = '#'

    for line in open(rocketfuel_path, "r").readlines():
        if comment_char in line:
            # split on comment char, keep only the part before
            line, _ = line.split(comment_char, 1)
            line = line.strip()
        if len(line) == 0:
            continue
        # Parse line.
        if line.startswith("-"):
            # Case external node
            # -euid =externaladdress rn
            try:
                node = int(re.findall("-\d+", line)[0])
                address = (re.findall("=\S+", line)[0])[1:]  # .strip("=")
                r = int(re.findall("r\d$", line)[0][1:])  # .strip("r"))
            except IndexError:
                raise ValueError('Invalid input file. Parsing failed '\
                                 'while trying to parse an external node')

            topology.add_node(node, type='external', address=address, r=r)

        else:
            # Case internal node
            # uid @loc [+] [bb] (num_neigh) [&ext] -> <nuid-1> <nuid-2>
            # ... {-euid} ... =name[!] rn
            try:
                node = int(re.findall("\d+", line)[0])
                node_location = re.findall("@\S*", line)[0]
                node_location = re.sub("[\+@]", "", node_location)
                r = int(re.findall("r\d$", line)[0][1:])  # .strip("r"))
                address = (re.findall("=\S+", line)[0])[1:]  # .strip("=")
            except IndexError:
                raise ValueError('Invalid input file. Parsing failed '\
                                 'while trying to parse an internal node')
            internal_links = re.findall("<(\d+)>", line)
            external_links = re.findall("{(-?\d+)}", line)
            backbone = True if len(re.findall("\sbb\s", line)) > 0 \
                       else False
            topology.add_node(node, type='internal',
                              location=node_location,
                              address=address, r=r, backbone=backbone)
            for link in internal_links:
                link = int(link)
                if node != link:
                    topology.add_edge(node, link, type='internal')
            for link in external_links:
                link = int(link)
                if node != link:
                    topology.add_edge(node, link, type='external')
    return topology

if __name__ == "__main__":

    # use an ArgumentParser for a nice CLI
    parser = argparse.ArgumentParser()

    # options (self-explanatory)
    parser.add_argument(
        "--mode", 
         help = """type of topologies to work with. e.g. 'pop' for pop-level 
            topologies, 'isp-map' for isp maps, etc.""")

    parser.add_argument(
        "--data-path", 
         help = """path w/ topology files (e.g. .cch for isp maps, edge folders 
            for pop-level maps, etc.)""")

    args = parser.parse_args()

    # quit if a dir w/ causality files hasn't been provided
    if not args.data_path:
        sys.stderr.write("""%s: [ERROR] please supply a data path\n""" % sys.argv[0]) 
        parser.print_help()
        sys.exit(1)

    # extract the contents of the .cch file to a networkx object
    if args.mode == "pop":
        topology = parse_pop_level_isp_map(args.data_path)
        convert_to_scn(topology)

    elif args.mode == "isp-map":
        topology = parse_isp_map(args.data_path)
        convert_to_scn(topology)

    else:
        sys.stderr.write("""%s: [ERROR] invalid mode (%s)\n""" % (sys.argv[0], args.mode)) 
        parser.print_help()
        sys.exit(1)

    # draw the graph
    nx.draw_spring(topology)
    # save the figure in <rocketfuel-file>.pdf
    plt.savefig(args.data_path.rstrip("/").rstrip(".cch").split("/")[-1].split(":")[-1] + ".pdf", bbox_inches='tight', format = 'pdf')


