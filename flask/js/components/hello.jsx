import React from 'react';

class Hello extends React.Component {
  render() {
    return(
      <div>
        <label>
          Upload MAF Alignment
          <input type="file" />
        </label>
        <label>
          Upload New Genome (FASTA)
          <input type="file" />
        </label>
        <label>
          Upload New Genome (FASTA)
          <input type="file" />
        </label>
        <label>
          Upload Reference Genome (FASTA)
          <input type="file" />
        </label>
        <label>
          N-length
          <input type="number" />
        </label>
        <label>
          Combine Method
          <select>
            <option>First</option>
            <option>Best</option>
            <option>Longest</option>
          </select>
        </label>
      </div>
    )
  }
}

export default Hello;