from flask import Flask, render_template, jsonify

app = Flask(__name__)

@app.route("/")

def index():
  return render_template('index.html')

# def genomeBuild():
#   # inputs: refernce genome (1 "fasta" file)
#   # outputs: A bunch of files

#   # options:

# def lastAlign():
#   # inputs:

def mafToBed():
  return({})
  # INPUT
  # alignment in maf format 1 file - maf
  # new genome - 1 file (fasta)

  # OUTPUT DIR
  # new file for each sequence (bed)
  # a file *allOverviews.bed

def overViewToIntersect():
  return({})
  # INPUT
  # *allOverviews.bed
  # map (.bed file) - can be created by ref genome, number: makeGenomeMap.makeInexactGenomeMap || makeGenomeMap.makeExactGenomeMap
  # Output folder from before

  # INPUT
  # bed dir

def reOrderScaffolds():
  return({})
  # INPUT
  # bed dir, map (.bed file, same as above), ref genome, new genome,
  # options:
  # name of extra pieces group (ungroupedChrom)
  # combine method ("first", "best", "longest")
  # consecutiveOnly (bool)

if __name__ == "__main__":
  app.run()