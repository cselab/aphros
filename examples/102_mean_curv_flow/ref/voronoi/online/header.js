function loadCSS(filename) {
	var fileref=document.createElement("link");
	fileref.setAttribute("rel", "stylesheet");
	fileref.setAttribute("type", "text/css");
	fileref.setAttribute("href", filename);
	document.getElementsByTagName("head")[0].appendChild(fileref);
}

function loadHeader() {
	loadCSS("http://fonts.googleapis.com/css?family=Molengo");
	loadCSS("http://fonts.googleapis.com/css?family=Cantarell");
	//loadCSS("http://fonts.googleapis.com/css?family=Inconsolata");
	loadCSS("http://alexbeutel.com/styles.css");
	//loadCSS("/styles.css");
	
	var div2 = document.createElement('div');
	div2.id = 'header2';
	document.body.appendChild(div2);

	var div = document.createElement('div');
	div.id = 'header';
	
	var span = document.createElement('span');
	span.id='name';
	span.innerHTML = "Alex Beutel<span style='letter-spacing:0.075em;padding-left: 3px;'>...</span>";
	div.appendChild(span);
	
	span = document.createElement('span');
	span.id='links';
	span.innerHTML = "<a href='http://alexbeutel.com'>About Me</a><a href='http://alexbeutel.com/projects.html'>Research</a><a style='padding-right:10px' href='http://alexbeutel.com/AlexBeutelCV.pdf'>CV</a>";
	div.appendChild(span);
	
	document.body.appendChild(div);
	//document.body.insertBefore(div, document.body.children[0]);
	
}
loadHeader();
