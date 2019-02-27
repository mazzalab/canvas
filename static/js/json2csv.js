function JSONToCSVConverter(JSONData, ReportTitle, ShowLabel, cnv, category, position_input) {
  var arrData = JSONData //typeof JSONData != 'object' ? JSON.parse(JSONData) : JSONData;
  var CSV = '';
//Set Report title in first row or line

  console.log(ShowLabel)
  console.log(cnv)
  console.log(position_input)
  console.log(category)
  console.log(JSONData)

  //This condition will generate the Label/Header
  if (ShowLabel) {
      var row = "";
      //This loop will extract the label from 1st index of on array
      for (var position in arrData[category]) {
        if (position == position_input){
          CSV += ReportTitle + ' for '+arrData[category][position][cnv]['cnv']+'\r\n\n';
          for (var index in arrData[category][position][cnv]['data'][0]) {
          //Now convert each value to string and comma-seprated
            row += index + ',';
          }
          row = row.slice(0, -1);
          CSV += row + '\r\n';
        }
      }
  }

  // If there's a single element in the JSON, cannot iterate with for.
  if(!arrData[category][position_input][cnv]['data'].length && arrData[category][position_input][cnv]['data']){
    var row = "";
    for (var index in arrData[category][position_input][cnv]['data']) {
      row += '"' + arrData[category][position_input][cnv]['data'][index] + '",';
    }
    row.slice(0, row.length - 1);
    CSV += row + '\r\n';
  }

  // If JSON contains more than one element, a normal for can be used
  for (var i = 0; i < arrData[category][position_input][cnv]['data'].length; i++) {
    var row = "";
    for (var index in arrData[category][position_input][cnv]['data'][i]) {
      row += '"' + arrData[category][position_input][cnv]['data'][i][index] + '",';
    }
    row.slice(0, row.length - 1);
    CSV += row + '\r\n';
  }

  if (CSV == '') {
    alert("Invalid data");
    return;
  }

  var fileName = "Report_";
  fileName += ReportTitle.replace(/ /g,"_");
  var uri = 'data:text/txt;charset=utf-8,' + escape(CSV);
  var link = document.createElement("a");
  link.href = uri;

  //set the visibility hidden so it will not effect on your web-layout
  link.style = "visibility:hidden";
  link.download = fileName + ".txt";

  //this part will append the anchor tag and remove it after automatic click
  document.body.appendChild(link);
  link.click();
  document.body.removeChild(link);
}
