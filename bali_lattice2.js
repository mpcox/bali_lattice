  var a = 0.5;
  var b = 9.6;
  var restarter = 0;
  var pauser = 1;
  var ender = 0;
  var t = 0;
  var harvest_record = [];
  var alpha = 0;
  
  // MPC
  // Now allows 'run/pause' button to toggle pauser between 0 and 1
  function run_pause() {
     if (pauser == 0) {
     	pauser = 1;
     } else {
     	pauser = 0;
     }
     // console.log("pauser button pressed");
     // console.log("pauser is " + pauser);
  }
  
  // MPC
  // Add a synchronous sleep function (not natively in javascript)
  // Synchronous = code stops until this sleep function completes
  // Must be called with functions with 'async' keyword
  function sleep(ms) {
      return new Promise(resolve => setTimeout(resolve, ms));
  }
  
  function main()
  {
    // Geometry of the simulation domain
    var cells_x = 100; // Number of Cells in x direction
    var cells_y = 100; // Number of cells in y direction

  // Geometry of the simulation domain
    var ns = 4
    var Rpest = 2
    var Rharvest = 1
    var noise = 0.05
    var noiseblock = 4
    var H0 = 20.0
    
    // Colors
    var col0 = "rgb(10,153,0)"; 
    var col1 = "rgb(255,255,51)";     
    var col2 = "rgb(0,0,255)";     
    var col3 = "rgb(255,0,0)"; 
    var colArray = [col0, col1, col2, col3];

    var spin = create_non_unique_random_array(cells_y,cells_x,1,(ns))
    var pest = create_zeros_array(cells_y,cells_x);
    var water = create_zeros_array(cells_y,cells_x);
    var harvest = create_zeros_array(cells_y,cells_x);
    var fq = [];

    var cv = document.getElementById("cvMain");
    var ctx = cv.getContext("2d");
    var cv2 = document.getElementById("cvSecond");
    var cv3 = document.getElementById("cvThird");

    init();

    //////////////////////////////////////////////////

    function init() {

      var spin = create_non_unique_random_array(cells_y,cells_x,1,(ns))

      plot_spin(spin)

      Run_Lattice_Model();
      
    }
   
    /////////////////////////////////////////////////////
    
    // MPC
    // New runtime logic, checks for pauser and restarter
    async function Run_Lattice_Model() {
      
      for (var it = 0; it < 100000; it++) {
                    
          // restarter
          if(restarter == 1){
              spin = create_non_unique_random_array(cells_y,cells_x,1,(ns));
              plot_spin(spin);
              restarter = 0;
          }
          
          if(pauser == 1){
              await sleep(100);
          } else {
              window.setTimeout(IterateLattice, 100);
              await sleep(500);
          }
          
          // console.log("iteration is " + it + " and pauser is " + pauser);
      }

    }
    
    // MPC
    // This was the old version
    
    // function Run_Lattice_Model() {
    // 
    //   for (var it = 0; it < 10000; it++) {
    //       window.setTimeout(IterateLattice,500*it);
    //   }
    // 
    // }

    /////////////////////////////////////////////////////
	
    function IterateLattice() {
      
      // MPC
      // pauser and restarter logic used to be here
      
      // if (restarter > 0) {    
      //   restarter = 0;
      //   spin = create_non_unique_random_array(cells_y,cells_x,1,(ns));
      // }

      // if (pauser > 0) {
      //   pauser = 0;
      //   alert("Paused.");
      // }
      
      var spin2 = create_zeros_array(cells_x,cells_y);
      var water = create_zeros_array(cells_x,cells_y);
      var pest = create_zeros_array(cells_x,cells_y);
      var harvest = create_zeros_array(cells_x,cells_y);

      for (var j=0; j<ns; j++) {
        fq[j] = array_value_counter(spin,(j+1))/(cells_x*cells_y);
      }

      for (var ii=0; ii<cells_x; ii++) {
        for (var jj=0; jj<cells_y; jj++) {
          var ilimit1 = Math.max(0, (ii - Rpest));
          var ilimit2 = Math.min((cells_x - 1), (ii + Rpest));
          var Sneigh = [];
          for (var xx = ilimit1; xx < (ilimit2 + 1); xx++) {
            var width = (Rpest - Math.abs(xx - ii));
            var jlimit1 = Math.max(0, (jj - width));
            var jlimit2 = Math.min((cells_y - 1), (jj + width));
            for (var yy = jlimit1; yy < (jlimit2+1); yy++) {
              Sneigh.push(spin[xx][yy]);    
            }
          }
          water[ii][jj] = fq[(spin[ii][jj] - 1)];
          var f =  1.0 / (0.1 + (vector_value_counter(Sneigh, spin[ii][jj]) - 1.0)/(Sneigh.length - 1.0));
          var f2 = ( 1.0 / (0.1 + (vector_value_counter((Sneigh,(spin[ii][jj] - 1.0)))/(Sneigh.length - 1.0)) ) )
          pest[ii][jj] = ( 1.0 / (0.1 + (vector_value_counter(Sneigh, spin[ii][jj]) - 1.0)/(Sneigh.length - 1.0)) );
          harvest[ii][jj] = H0 - (a*pest[ii][jj]) - (b*water[ii][jj]);
        }
      }

      harvest_record.push(array_sum(harvest)/Math.pow((cells_x+0.0),2));

      for (var ww = 0; ww < cells_x; ww++) {
        for (var zz = 0; zz < cells_y; zz++) {
          var wlimit1 = Math.max(0, (ww-Rharvest));
          var wlimit2 = Math.min((cells_x - 1), (ww + Rharvest));
          var Hneigh = [];
          var SSneigh = [];
          for (var qx = wlimit1; qx < (wlimit2 + 1); qx++) {
            var width = Rharvest - Math.abs(qx - ww);
            var zlimit1 = Math.max(0, (zz - width));
            var zlimit2 = Math.min((cells_y - 1), (zz + width));
            for (qy = zlimit1; qy < (zlimit2 + 1); qy++) {
              SSneigh.push(spin[qx][qy]);
              Hneigh.push(harvest[qx][qy]);
            }
          }
          mh = Math.max.apply(null, Hneigh);
          Hmaxind = [];
          for (ms = 0; ms < (Hneigh.length); ms++) {
            if (Hneigh[ms] == mh) {
              Hmaxind.push(ms);
            }
          }
          var Hmi = Hmaxind[Math.floor(Math.random()*Hmaxind.length)];
          if (Hneigh[Hmi] > harvest[ww][zz]) {
            spin2[ww][zz] = SSneigh[Hmi];
          } else {
            spin2[ww][zz] = spin[ww][zz];
          }
        }
      }

      if (noise > 0) {
        var noise_counter = 0;
        while (noise_counter < noise*(cells_x*cells_y)) {
          inb = Math.floor(Math.random()*cells_x);
          jnb = Math.floor(Math.random()*cells_y);
          snb = 1 + Math.floor(Math.random()*(noiseblock));
          if (((inb + snb) <= cells_x)&&((jnb+snb) <= cells_y)) {
            var nss = ( 1 + Math.floor(Math.random()*(ns)) );
            for (var sx=0; sx<snb; sx++) {
              for (var sy=0; sy<snb; sy++) {
                spin2[inb+sx][jnb+sy] = nss;
              }
            }
          }
          noise_counter = noise_counter + Math.pow(snb,2);
        }
      }

      spin = spin2;

      var PatchSizeOut = PatchSize();
      S = PatchSizeOut[0];
      T = PatchSizeOut[1];
      cluster = PatchSizeOut[2];
      SlopeOut = Slope();

      plot_spin(spin)
      plot_harvest(harvest_record)
      plot_patchsize(SlopeOut)
      
      // MPC
      // The following lines returns NULL
      // As far as I can tell, it does nothing, as there is no 'demo' in the html
      
      // document.getElementById("demo").innerHTML = t;

      // t = t + 1;
    }

 /////////////////////////////////////////////////////

    function random_integer(min,max) {
      return (Math.round((max-min)*Math.random() + min));
    }

 /////////////////////////////////////////////////////

    function create_non_unique_random_array(num_rows,num_cols,min,max) {

    var matrix = [];

    for (var rowindex=0; rowindex<num_rows; rowindex++) {

      var row = [];

      for (var element=0; element<num_cols; element++) {
        row[element] = random_integer(min,max);
      }
      matrix.push(row) 
    }
    return (matrix);
    }

  /////////////////////////////////////////////////////   

    function create_zeros_array(num_rows,num_cols) {

    var matrix = [];

    for (var rowindex=0; rowindex<num_rows; rowindex++) {

      var row = [];

      for (var element=0; element<num_cols; element++) {
        row[element] = 0;
      }
      matrix.push(row) 
    }
    return (matrix);
    }

 /////////////////////////////////////////////////////

    function create_negones_vector(vlength) {

    var vector = [];
    for (var element=0; element<vlength; element++) {
      vector[element] = -1;
    } 
    return (vector);
    }

 /////////////////////////////////////////////////////

    function create_zeros_vector(vlength) {

    var vector = [];
    for (var element=0; element<vlength; element++) {
      vector[element] = 0;
    } 
    return (vector);
    }

 /////////////////////////////////////////////////////

    function create_negones_array(num_rows,num_cols) {

    var matrix = [];
    for (var rowindex=0; rowindex<num_rows; rowindex++) {
      var row = [];
      for (var element=0; element<num_cols; element++) {
        row[element] = -1;
      }
      matrix.push(row)
    } 
    return (matrix);
    }

 /////////////////////////////////////////////////////

    function array_value_counter(array,value) {

    var counter = 0;

    for (var rowindex=0; rowindex<array.length; rowindex++) {
      for (var element=0; element<array[rowindex].length; element++) {
        if (array[rowindex][element] == value) {
          counter = counter + 1;
        }
      }
    }
    return (counter);
    }

 /////////////////////////////////////////////////////

    function vector_value_counter(array,value) {

    var counter = 0;
    for (var element=0; element<array.length; element++) {
      if (array[element] == value) {
        counter = counter + 1;
      }
    }
    return counter
    }

 /////////////////////////////////////////////////////

    function vector_notvalue_counter(array,value) {

    var counter = 0;
    //for (var rowindex=0; rowindex<array.length; rowindex++) {
    for (var element=0; element<array.length; element++) {
      if (array[element] != value) {
        counter = counter + 1;
      }
    }
    return (counter);
    }

 /////////////////////////////////////////////////////

    function array_sum(array) {

    var sum = 0;

    for (var rowindex=0; rowindex<array.length; rowindex++) {
      for (var element=0; element<array[rowindex].length; element++) {
        sum = sum + array[rowindex][element];
      }
    }
    return (sum);
    }

 /////////////////////////////////////////////////////

    function vector_sum(vector) {

    var sum = 0;

    for (var element=0; element<vector.length; element++) {
      sum += vector[element];
    }
    return sum;
    }

 /////////////////////////////////////////////////////

    function plot_spin(array) {
    for (var ii=0; ii<array.length; ii++) {
      for (var jj=0; jj<array[ii].length; jj++) {
        ctx.fillStyle = colArray[spin[ii][jj]-1];
        ctx.fillRect(ii*(cv.width/cells_x), jj*(cv.height/cells_y), (cv.width/cells_x), (cv.height/cells_y));
      }
    }
    }

 /////////////////////////////////////////////////////

    function plot_harvest(vector) {
      var plotObject = {
      x: time_axis(vector.length),
      y: vector,
      mode: 'lines',
      line: {
        width: 1
      } };

      Plotly.newPlot( cv2, [plotObject], {
      width: 400,
      height: 230,
      title: 'Harvest',
      xaxis: {
        title: 'Iteration',
        titlefont: {
          size: 12,
        }
      },
      yaxis: {
        title: 'Harvest',
        titlefont: {
           size: 12,
        },
        tickfont:{size:8}
      },
  
      margin: { l:34, r:34, b: 34, t: 34, pad: 4 } }, {displayModeBar: false} );

    }

 /////////////////////////////////////////////////////

    function time_axis(time_current) {
      var times = [];
      for (var i = 0; i < time_current; i++) {
        times.push(i);
      }
      return times;
    }

 /////////////////////////////////////////////////////

    function PatchSize() {
      var N = cells_y;
      var M = cells_x;
      var k = 4;
      var Alist = create_negones_array(N*M,k);
      for (var jj = 0; jj < N; jj++) {
        for (var kk = 0; kk < M; kk++) {
          var v = (kk*N) + jj;
          var link = 0;
          if ((jj - 1) > -1) {
            var vN = (kk*N) + (jj - 1);
            Alist[v][link] = vN;
            link += 1;
          }
          if ((jj + 1) < N) {
            var vS = (kk*N) + (jj + 1);
            Alist[v][link] = vS;
            link += 1;
          }
          if ((kk + 1) < M) {
            var vE = ((kk + 1)*N) + jj;
            Alist[v][link] = vE;
            link += 1;
          }
          if ((kk - 1) < -1) {
            var vW = ((kk - 1)*N) + jj;
            Alist[v][link] = vW;
            link += 1;
          }
        }
      }

      var degree = create_zeros_vector(N*M);

      for (var ii = 0; ii < (N*M); ii++) {
        degree[ii] = vector_notvalue_counter(Alist[ii], -1);
      }

      var cluster = create_negones_vector((N*M));
      var Nc = 0;
      var T = create_zeros_vector(N*M);

      for (var nn = 0; nn < (N*M); nn++) {

        var origin = nn;

        if (cluster[origin]==(-1)) {

          var crop = spin[(Math.floor(origin/N))][origin % N];

          if (crop > 0) {

            cluster[origin] = Nc;

            T[Nc] = crop;
            var queue = create_negones_vector(N*M);
            queue[0] = origin;
            var q = 0;
            var uu = 0;
            while ((queue[uu] > (-1))&&(uu<(N*M))) {
              var current = queue[uu];
              //alert("queue is " + queue)
              var connected = []
              for (var rr = 0; rr < degree[current]; rr++) {
                connected.push(Alist[current][rr]);
              }

              for (var j = 0; j < connected.length; j++) {

                if ((spin[Math.floor(connected[j]/N)][connected[j] % N] == crop)&&(cluster[connected[j]] == (-1))) {
                  cluster[connected[j]] = Nc;
                  queue[(q + 1)] = connected[j];
                  q += 1;
                }
              }
              uu += 1;
            }
            Nc += 1;
          }    
        }
        }

    var T0 = [];
    for (var tt = 0; tt <= (Math.max.apply(null, cluster) + 1) ; tt++) {
      T0.push(T[tt]);
    }
    var S = create_zeros_vector((Math.max.apply(null, cluster) + 1));
    for (var cc = 0; cc < (Math.max.apply(null, cluster) + 1); cc++) {
      S[cc] = vector_value_counter(cluster, cc);
    }
    return [S, T, cluster]
    }

 /////////////////////////////////////////////////////

    function Slope() {

      var c = 3.0;
      var x_min = Math.pow( (Math.sqrt(3200.0)/20.0), 2);
      var power1 = Math.log(x_min)/Math.log(c);
      var x_fit = vectorpower(c, arangea2(power1, Math.ceil(Math.log(cells_x*cells_y)/Math.log(c))));
      var y_fit0 = [];
     // alert('xmin ' + x_min)
      ///alert('power1 ' + power1)
      //alert('xfit ' + x_fit)

      for (var szi = 0; szi < x_fit.length; szi++) {
        var sz = x_fit[szi];
        y_fit0.push(greaterthancounter(S,sz))
      }

      //alert('y_fit0 ' + y_fit0)

      var y_fit = []
      for (var ind = 0; ind < y_fit0.length; ind++) {
        y_fit[ind] = y_fit0[ind]/y_fit0[0];
      }

      //alert('y_fit ' + y_fit)

      var bmin = 1;
      var LL = greaterthancounter(S, x_fit[bmin]);

      //alert(LL)

      // var A1 = dotprod(truncatevector(y_fit,bmin-1),truncatevector(arangea(1,y_fit.length+1),bmin-1));
      // var A2 = ((power1 - 1) - (Math.log(x_fit[bmin - 1])/Math.log(c)) + A1 );
      // var A3 = (1/A2);
      // var A4 = Math.log(1 + A3);
      // var alpha_hat = 1 + ( A4 / Math.log(c));

      var A4 = dotprod(y_fit,(arangea(1,y_fit.length+1))) / vector_sum(y_fit);
      var A3 = (power1 - 1.0) - (Math.log(x_fit[0]) / Math.log(c)) + A4;
      var A2 = Math.log(1 + (1.0 / A3));
      var A1 = A2 / Math.log(c);
      var alpha_hat = 1 + A1;
      
      
      
      

      
      var sigma = ((Math.pow(c,alpha_hat))-c)/(Math.pow(c,(0.5*(1+alpha_hat)))*Math.log(c)*Math.sqrt(LL));
      var alpha = alpha_hat - 1;

      return [x_fit, y_fit, alpha, sigma, bmin, c]
    }

 /////////////////////////////////////////////////////

    function dotprod(vec1,vec2) {
      var dot = 0;
      for (var element=0; element < vec1.length; element++) {
        dot += vec1[element]*vec2[element];
      }
      return dot
    }

 /////////////////////////////////////////////////////

    function truncatevector(vector,minindex) {
      var vectort = [];
      counter = 0;
      for (var vi=minindex; vi < vector.length; vi++) {
        vectort[counter] = vector[vi];
        counter += 1;
      }
      return vectort
    }

 /////////////////////////////////////////////////////

    function greaterthanvector(vector,value) {
      var vectorgt = [];
      for (var vi=0; vi < vector.length; vi++) {
        if (vector[vi] >= value) {
          vectorgt.push(vector[vi]);
        }
      }
      return vectorgt
    }

 /////////////////////////////////////////////////////

    function greaterthancounter(vector,value) {
      var count = 0;
      for (var vi=0; vi < vector.length; vi++) {
        if (vector[vi] >= value) {
          count += 1;
        }
      }
      return count
    }

 /////////////////////////////////////////////////////

    function arangea(left,right) {
      var vector = [];
      var element = Math.ceil(left);
      var counter = 0;
      while (element < right) {
        vector[counter] = element;
        counter += 1;
        element += 1;
      }
      return vector
    }

 /////////////////////////////////////////////////////

    function arangea2(left,right) {
      var vector = [];
      var element = (left);
      var counter = 0;
      while (element < right) {
        vector[counter] = element;
        counter += 1;
        element += 1;
      }
      return vector
    }

 /////////////////////////////////////////////////////

    function vectorpower(value, vector) {
      var vpow = [];
      for (var element = 0; element < vector.length; element++) {
        vpow[element] = Math.round(Math.pow(value, vector[element]));
      }
      return vpow
    }

 /////////////////////////////////////////////////////

    function plot_patchsize(SlopeOut) {
      x_fit = SlopeOut[0];
      y_fit = SlopeOut[1];
      alpha = SlopeOut[2];
      sigma = SlopeOut[3];
      bmin = SlopeOut[4];
      c = SlopeOut[5];

      var x0 = arangea(Math.pow(10,0.7), Math.pow(10,3.1));
      var y0 = [];

      for (var ind=0; ind < x0.length; ind++) {
        y0[ind] = Math.pow(x0[ind], (-alpha))/Math.pow(x_fit[0], (-alpha));
      }

      var plotObject = {
        x: x_fit,
        y: y_fit,
        name: 'Patch Size Distribution',
        mode: 'markers',
        };

      var plotObject2 = {
        x: x0,
        y: y0,
        name: 'Cumulative Distribution Function',
        mode: 'lines',
      }

        Plotly.newPlot( cv3, [plotObject, plotObject2], {
        title: 'Patch Size Distribution',
        width: 400,
        height: 230,
        margin: 50,
        xaxis: {
          type: 'log',
          range: [0.6, 3.1],
          title: 'Size',
        },
        yaxis: {
          type: 'log',
          range: [-2.3, 0.4],
          title: 'cumulative distribution function',
          titlefont:{size:8},
          tickfont: {size:8},
        },
        showlegend: false,

        annotations: [ {x: 1.2, y: -1.5, text: '\u03B1 = ' + alpha.toString().slice(0,6), font: {size:20}, showarrow: false}],
        margin: { l:34, r:34, b: 34, t: 34, pad: 4 } }, {displayModeBar: false});
      }



}
