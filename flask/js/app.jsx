import React from 'react';
import ReactDOM from 'react-dom';

import Hello from './components/hello';

document.addEventListener('DOMContentLoaded', () => {
  const root = document.getElementById('reactEntry');
  ReactDOM.render(<Hello/>, root);
})