var chart = document.getElementById("chart");
var margin = {
  bottom: 80,
  top: 20,
  left: 50,
  right: 30,
};
var width = chart.clientWidth - margin.left - margin.right;
var height = 400 - margin.top - margin.bottom;

var svg = d3.select("#chart")
  .append("svg")
    .attr("width", width + margin.left + margin.right)
    .attr("height", height + margin.top + margin.bottom)
  .append("g")
    .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

var data = gData;

rev2row = {};
for (d of data) {
  rev2row[d.rev] = d;
}

var barshift = width / data.length
var barwidth = barshift * 0.8

var x = d3.scaleBand()
    .domain(data.map(x => x.rev))
    .range([0, width])
svg.append("g")
    .attr("transform", "translate(0," + height + ")")
    .call(d3.axisBottom(x))
    .selectAll("text")
      .attr("transform", "translate(10,5)rotate(-45)")
      .attr("class", "xticklabel")
      .text(function(d) { return rev2row[d].date; })

var y = d3.scaleLinear()
    .range([height, 0]);
    y.domain([
      d3.min(data, function(d) { return d.runtime; }) * 0.8,
      d3.max(data, function(d) { return d.runtime; }) * 1.05
    ]);
svg.append("g")
    .call(d3.axisLeft(y))
    .selectAll("text")
      .attr("class", "yticklabel")
svg.append("text")
    .attr("transform", "rotate(-90)")
    .attr("y", -margin.left)
    .attr("x", -(height / 2))
    .attr("dy", "1em")
    .attr("class", "ylabel")
    .text("runtime [s]");

var tooltip = d3.select("body").append("div")
    .attr("class", "tooltip")
    .style("opacity", 0);

svg.selectAll("rect")
    .data(data)
    .enter()
    .append("svg:a").attr("xlink:href", function(d) { return d.url; })
    .append("rect")
      .attr("class", "bar")
      .attr("transform", function(d) {
        xpos = x(d.rev) + (barshift - barwidth) * 0.5;
        ypos = y(d.runtime);
        return "translate(" + xpos + "," + ypos + ")";
      })
      .attr("width", barwidth)
      .attr("height", function(d) { return height - y(d.runtime); })

svg.selectAll("rect")
    .on("mouseover", function(d){
      tooltip.transition()
          .duration(50)
          .style("opacity", .9);
      t = parseFloat(d.runtime).toFixed(2)
      tooltip.html(t + " s<br/>"  + d.date + "<br/>"  + d.rev)
          .style("left", (d3.event.pageX) + "px")
          .style("top", (d3.event.pageY - 50) + "px");
    })
    .on("mousemove",function(){
      tooltip
        .style("left", (d3.event.pageX) + "px")
        .style("top", (d3.event.pageY - 50) + "px");
    })
    .on("mouseout",function(){
      tooltip.transition()
          .duration(50)
          .style("opacity", 0);
    })

