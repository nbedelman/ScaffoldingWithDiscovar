var webpack = require('webpack');
module.exports = {
  devtool: "source-map",
  context: __dirname,
  entry: [
    "./js/app.jsx"
  ],
  output: {
    path: __dirname + '/static',
    filename: 'bundle.js'
  },
  module: {
    loaders: [
      {
        test: /\.(js|jsx)$/,
        loader: 'babel-loader',
        query: {
          presets: ['es2015','stage-2','react']
        },
        exclude: /node_modules/
      }
    ]
  },
  resolve: {
    extensions: ['.js', '.jsx'],
  }
};