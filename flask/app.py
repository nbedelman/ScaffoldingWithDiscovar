from flask import Flask, render_template, jsonify

app = Flask(__name__)

@app.route("/")

def index():
    return render_template('index.html')

@app.route("/api/v1/get_word", methods=['GET'])

def index():
    return jsonify({"foo":"bar")

if __name__ == "__main__":
    app.run()