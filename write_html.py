import json
names = ["B2Mg1-108064"]

json_file = open("data.json", "r")
data_json = json.load(json_file)

for name in data_json.keys():

    eform = data_json[name]["eform"]
    dosf = data_json[name]["dosf"]
    magt = data_json[name]["magt"]
    atom_type = data_json[name]["type"]
    nat = 0
    for atom in atom_type.keys():
        nat += int(atom_type[atom]["nat"])

    with open("/home/kawamura/work/xsf_icsd/" + str(name) + ".xsf") as f:
        lattice = f.readlines()
        lattice = "".join(lattice).replace("\n", "\\n")

    with open("html/" + name + ".html", 'w') as f:
        print("<html>", file=f)
        print("  <head>", file=f)
        print("    <title>" + name + "</title>", file=f)
        print("    <meta charset=\"utf-8\">", file=f)
        print("    <script type=\"text/javascript\" src=\"../JSmol.min.js\"></script>", file=f)
        print("    <script type=\"text/javascript\"> ", file=f)
        print("// supersimple2.htm - illustrating the use of jQuery(document).ready to ", file=f)
        print("// populate all spans and divs AFTER the page is loaded.", file=f)
        print("// This is good programming practice.", file=f)
        print("$(document).ready(", file=f)
        print("function() {", file=f)
        print("Info = {", file=f)
        print("	width: 600,", file=f)
        print("	height: 600,", file=f)
        print("	debug: false,", file=f)
        print("	j2sPath: \"../j2s\",", file=f)
        print("	color: \"0xC0C0C0\",", file=f)
        print("  disableJ2SLoadMonitor: true,", file=f)
        print("  disableInitialConsole: true,", file=f)
        print("	addSelectionOptions: false,", file=f)
        print("	serverURL: \"https://chemapps.stolaf.edu/jmol/jsmol/php/jsmol.php\",", file=f)
        print("	use: \"HTML5\",", file=f)
        print("	readyFunction: null,", file=f)
        print("    script: \"load INLINE \\\"" + lattice +
              "\\\" {3,3,3}; set perspectiveDepth ON; select;set defaultLabelXYZ \\\"%e\\\"; set labelToggle; "
              + "set labelAtom; set labelOffset 0 0; color labels black\"", file=f)
        print("}", file=f)
        print("$(\"#mydiv\").html(Jmol.getAppletHtml(\"jmolApplet0\",Info))", file=f)
        print("}", file=f)
        print(");", file=f)
        print("    </script>", file=f)
        print("<!-- gnuplot head start -->", file=f)
        print("    <meta http-equiv=\"content-type\" content=\"text/html; charset=UTF-8\">", file=f)
        print("    <!--[if IE]><script type=\"text/javascript\" src=\"excanvas.js\"></script><![endif]-->", file=f)
        print("    <script src=\"../gnuplot/canvastext.js\"></script>", file=f)
        print("    <script src=\"../gnuplot/gnuplot_common.js\"></script>", file=f)
        print("    <script src=\"../gnuplot/gnuplot_dashedlines.js\"></script>", file=f)
        print("    <script src=\"../gnuplot/gnuplot_mouse.js\"></script>", file=f)
        print("    <script type=\"text/javascript\"> gnuplot.help_URL = \"/canvas_help.html\"; </script>", file=f)
        print("    <script type=\"text/javascript\">", file=f)
        print("var canvas, ctx;", file=f)
        print("gnuplot.grid_lines = true;", file=f)
        print("gnuplot.zoomed = false;", file=f)
        print("gnuplot.active_plot_name = \"gnuplot_canvas\";", file=f)
        print("    </script>", file=f)
        print("    <script src=\"../band/" + name + ".js\"></script>", file=f)
        print("    <link type=\"text/css\" href=\"../gnuplot/gnuplot_mouse.css\" rel=\"stylesheet\">", file=f)
        print("	<!-- gnuplot head end -->", file=f)
        print("  </head>", file=f)
        print("  <body bgcolor=\"CCFFCC\" onload=\"gnuplot_canvas(); " +
              "gnuplot.init();\" oncontextmenu=\"return false;\">", file=f)
        print("	<h1>" + name + "</h1>", file=f)
        print("	<h2>Crystal structure</h2>", file=f)
        print("	<span id=mydiv></span>", file=f)
        print("  <h2>Properties</h2>", file=f)
        print("  <p>Number of atoms per unit cell : " + str(nat) + "</p>", file=f)
        print("  <p>Formation energy : " + str(eform) + " eV/atom</p>", file=f)
        print("  <p>Magnetization : " + str(magt) + " &mu;<sub>B</sub>/atom</p>", file=f)
        print("  <p>DOS at E<sub>F</sub> : " + str(dosf) + " /eV/atom (both spin)</p>", file=f)
        print("	<h2>Band structure</h2>", file=f)
        print("    <div class=\"gnuplot\">", file=f)
        print("      <canvas id=\"Tile\" width=\"32\" height=\"32\" hidden></canvas>", file=f)
        print("      <table class=\"mbleft\">", file=f)
        print("        <tr>", file=f)
        print("          <td class=\"mousebox\">", file=f)
        print("            <table class=\"mousebox\" border=0>", file=f)
        print("              <tr>", file=f)
        print("                <td class=\"mousebox\">", file=f)
        print("                  <table class=\"mousebox\" id=\"gnuplot_mousebox\" border=0>", file=f)
        print("                    <tr><td class=\"mbh\"></td></tr>", file=f)
        print("                    <tr>", file=f)
        print("                      <td class=\"mbh\">", file=f)
        print("                        <table class=\"mousebox\">", file=f)
        print("                          <tr>", file=f)
        print("                            <td class=\"icon\"></td>", file=f)
        print("                            <td class=\"icon\" onclick=gnuplot.toggle_grid>" +
              "<img src=\"../gnuplot/grid.png\" id=\"gnuplot_grid_icon\" class=\"icon-image\" " +
              "alt=\"#\" title=\"toggle grid\"></td>", file=f)
        print("                            <td class=\"icon\" onclick=gnuplot.unzoom>" +
              "<img src=\"../gnuplot/previouszoom.png\" id=\"gnuplot_unzoom_icon\" class=\"icon-image\" " +
              "alt=\"unzoom\" title=\"unzoom\"></td>", file=f)
        print("                            <td class=\"icon\" onclick=gnuplot.rezoom>" +
              "<img src=\"../gnuplot/nextzoom.png\" id=\"gnuplot_rezoom_icon\" class=\"icon-image\" " +
              "alt=\"rezoom\" title=\"rezoom\"></td>", file=f)
        print("                            <td class=\"icon\" onclick=gnuplot.toggle_zoom_text>" +
              "<img src=\"../gnuplot/textzoom.png\" id=\"gnuplot_textzoom_icon\" class=\"icon-image\" " +
              "alt=\"zoom text\" title=\"zoom text with plot\"></td>", file=f)
        print("                            <td class=\"icon\" onclick=gnuplot.popup_help()>" +
              "<img src=\"../gnuplot/help.png\" id=\"gnuplot_help_icon\" class=\"icon-image\" " +
              "alt=\"?\" title=\"help\"></td>", file=f)
        print("                          </tr>", file=f)
        print("                          <tr>", file=f)
        print("                            <td class=\"icon\" onclick=gnuplot.toggle_plot(\"gp_plot_1\")>1</td>",
              file=f)
        print("                            <td class=\"icon\" > </td>", file=f)
        print("                            <td class=\"icon\" > </td>", file=f)
        print("                            <td class=\"icon\" > </td>", file=f)
        print("                            <td class=\"icon\" > </td>", file=f)
        print("                            <td class=\"icon\" > </td>", file=f)
        print("                          </tr>", file=f)
        print("                        </table>", file=f)
        print("                      </td>", file=f)
        print("                    </tr>", file=f)
        print("                  </table>", file=f)
        print("                </td>", file=f)
        print("              </tr>", file=f)
        print("              <tr>", file=f)
        print("                <td class=\"mousebox\">", file=f)
        print("                  <table class=\"mousebox\" id=\"gnuplot_mousebox\" border=1>", file=f)
        print("                    <tr>", file=f)
        print("                      <td class=\"mb0\">x&nbsp;</td>", file=f)
        print("                      <td class=\"mb1\"><span id=\"gnuplot_canvas_x\">&nbsp;</span></td>", file=f)
        print("                    </tr>", file=f)
        print("                    <tr>", file=f)
        print("                      <td class=\"mb0\">y&nbsp;</td>", file=f)
        print("                      <td class=\"mb1\"><span id=\"gnuplot_canvas_y\">&nbsp;</span></td>", file=f)
        print("                    </tr>", file=f)
        print("                  </table>", file=f)
        print("                </td>", file=f)
        print("              </tr>", file=f)
        print("            </table>", file=f)
        print("          </td>", file=f)
        print("          <td>", file=f)
        print("            <table class=\"plot\">", file=f)
        print("              <tr>", file=f)
        print("                <td>", file=f)
        print("                  <canvas id=\"gnuplot_canvas\" width=\"600\" height=\"400\" tabindex=\"0\">", file=f)
        print("                    Sorry, your browser seems not to support the HTML 5 canvas element", file=f)
        print("                  </canvas>", file=f)
        print("                </td>", file=f)
        print("              </tr>", file=f)
        print("            </table>", file=f)
        print("          </td>", file=f)
        print("        </tr>", file=f)
        print("      </table>", file=f)
        print("    </div>", file=f)
        print("    <h2>Fermi surface and PDOS (/eV/atom both spin)</h2>", file=f)
        print("    <ul>", file=f)
        for atom in atom_type.keys():
            for n_main in atom_type[atom].keys():
                if not isinstance(atom_type[atom][n_main], dict):
                    continue
                for l_ang in atom_type[atom][n_main].keys():
                    orb = str(atom) + str(n_main) + str(l_ang)
                    pdos0 = sum(atom_type[atom][n_main][l_ang]["pdos"]) / nat
                    print("      <li>" + orb + "(<a href=\"../fermisurfer/index.php?frmsf=../fermi/"
                          + name + "_" + orb + "_1.js\" " + "target=\"_blank\">u</a>, " +
                          "<a href=\"../fermisurfer/index.php?frmsf=../fermi/"
                          + name + "_" + orb + "_2.js\" " + "target=\"_blank\">d</a>)"
                          + " : " + str(pdos0) + "</li>", file=f)
        print("    </ul>", file=f)
        print("    <h2>Input and output</h2>", file=f)
        print("    <ul>", file=f)
        print("      <li><a href=\"../scfin/scf_" + name + ".in\">SCF input</a></li>", file=f)
        print("      <li><a href=\"../scfout/scf_" + name + ".out\">SCF output</a></li>", file=f)
        print("    </ul>", file=f)
        print("  </body>", file=f)
        print("</html>", file=f)
