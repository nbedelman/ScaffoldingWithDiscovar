"""
SWD
Usage:
  SWD hello <d> <m> <fo> <fn>
  SWD lastAlign
  SWD mafToSwbed
  SWD reOrderScaffolds
  SWD -h | --help
  SWD --version
Options:
  <d>                               Directory with contig bed files
  <m>                               Original genome map file (AGP)
  <fo>                              Fasta file of reference genome
  <fn>                              Fasta file of new assembly
  -h --help                         Show this screen.
  --version                         Show version.
Help:
  For help using this tool, please open an issue on the Github repository:
  https://github.com/nbedelman/ScaffoldingWithDiscovar
"""


from inspect import getmembers, isclass

from docopt import docopt

from . import __version__ as VERSION


def main():
    """Main CLI entrypoint."""
    import commands
    options = docopt(__doc__, version=VERSION)

    # Here we'll try to dynamically match the command the user is trying to run
    # with a pre-defined command class we've already created.
    for k, v in options.iteritems():
        if hasattr(commands, k):
            module = getattr(commands, k)
            commands = getmembers(module, isclass)
            command = [command[1] for command in commands if command[0] != 'Base'][0]
            command = command(options)
            command.run()

