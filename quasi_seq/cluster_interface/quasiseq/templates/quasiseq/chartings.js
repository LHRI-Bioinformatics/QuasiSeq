var pieData = [
          {
              value: 20,
              color:"#878BB6"
          },
          {
              value : 40,
              color : "#4ACAB4"
          },
          {
              value : 10,
              color : "#FF8153"
          },
          {
              value : 30,
              color : "#FFEA88"
          }
      ];
      // Get the context of the canvas element we want to select
var countries= document.getElementById("myChart").getContext("2d");
new Chart(countries).Pie(pieData);
