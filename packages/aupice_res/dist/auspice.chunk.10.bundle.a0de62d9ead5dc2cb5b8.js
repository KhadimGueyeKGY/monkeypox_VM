(window.webpackJsonp=window.webpackJsonp||[]).push([[10],{1175:function(t,e,n){"use strict";n.r(e);var a=n(1),s=n.n(a),o=n(234),r=n.n(o),i=n(796),c=n(1157),l=n(817),h=(n(841),n(843)),p=n(865),u=n(822),d=n(946),f=n(842),y=n(923),v=n(1170),m=n(1148),x=n(1149),b=n(880),g=n.n(b),M=n(828),w=n(874),N=n(818);function O(t,e){var n="undefined"!=typeof Symbol&&t[Symbol.iterator]||t["@@iterator"];if(!n){if(Array.isArray(t)||(n=function(t,e){if(!t)return;if("string"==typeof t)return j(t,e);var n=Object.prototype.toString.call(t).slice(8,-1);"Object"===n&&t.constructor&&(n=t.constructor.name);if("Map"===n||"Set"===n)return Array.from(t);if("Arguments"===n||/^(?:Ui|I)nt(?:8|16|32)(?:Clamped)?Array$/.test(n))return j(t,e)}(t))||e&&t&&"number"==typeof t.length){n&&(t=n);var a=0,s=function(){};return{s:s,n:function(){return a>=t.length?{done:!0}:{done:!1,value:t[a++]}},e:function(t){throw t},f:s}}throw new TypeError("Invalid attempt to iterate non-iterable instance.\nIn order to be iterable, non-array objects must have a [Symbol.iterator]() method.")}var o,r=!0,i=!1;return{s:function(){n=n.call(t)},n:function(){var t=n.next();return r=t.done,t},e:function(t){i=!0,o=t},f:function(){try{r||null==n.return||n.return()}finally{if(i)throw o}}}}function j(t,e){(null==e||e>t.length)&&(e=t.length);for(var n=0,a=new Array(e);n<e;n++)a[n]=t[n];return a}var k=function(t,e,n,a,s){this.svg=Object(c.a)(t),this.annotations=e,this.geneMap=n,this.maxNt=a,this.callbacks=s,this.okToDrawBars=!1};k.prototype.render=function(t){this.props=t,this.aa="aa"===t.mutType,this.bars=t.bars,this.selectedNodes=Object(M.f)(t.colorBy)?this._getSelectedNodes(Object(M.a)(t.colorBy,t.geneLength)):[],this.svg.selectAll("*").remove(),this._calcOffsets(t.width,t.height),this._drawMainNavElements(),this._addZoomLayers(),this._setScales(this.maxNt+1,t.maxYVal),this.zoomCoordinates=Object(M.f)(t.colorBy)?this._getZoomCoordinates(Object(M.a)(t.colorBy,t.geneLength),t.geneMap):this.scales.xNav.domain(),this.zoomCoordinates[0]=t.zoomMin?t.zoomMin:this.zoomCoordinates[0],this.zoomCoordinates[1]=t.zoomMax?t.zoomMax:this.zoomCoordinates[1],this._drawAxes(),this._addBrush(),this._addClipMask(),this._drawGenes(this.annotations),this._drawZoomGenes(this.annotations),this.okToDrawBars=!0,this._drawBars(),this.zoomed=this._createZoomFn()},k.prototype.update=function(t){var e=this,n=t.aa,a=void 0===n?void 0:n,s=t.selected,o=void 0===s?void 0:s,r=t.newBars,i=void 0===r?void 0:r,c=t.maxYVal,l=void 0===c?void 0:c,h=t.clearSelected,p=void 0!==h&&h,u=t.gene,d=void 0===u?void 0:u,f=t.start,y=void 0===f?void 0:f,v=t.end,m=void 0===v?void 0:v,x=t.zoomMax,b=void 0===x?void 0:x,g=t.zoomMin,M=void 0===g?void 0:g,w=void 0!==a&&a!==this.aa;if((i||w)&&(w&&(this.aa=a),i&&(this.bars=i),this._updateYScaleAndAxis(l),this._drawBars()),void 0!==o?(this._clearSelectedBars(),this.selectedNodes=this._getSelectedNodes(o),this._highlightSelectedBars()):p&&this._clearSelectedBars(),void 0!==d&&void 0!==y&&void 0!==m){var O=m-y,j=d===N.K?0:1*O;this.navGraph.select(".brush").call(this.brush.move,(function(){return[Math.max(e.scales.xNav(y-j),e.scales.xNav(0)),Math.min(e.scales.xNav(m+j),e.scales.xNav(e.scales.xNav.domain()[1]))]}))}if(void 0!==M||void 0!==b){var k=void 0===M?0:M,C=void 0===b?this.scales.xNav.domain()[1]:b;this.navGraph.select(".brush").call(this.brush.move,(function(){return[e.scales.xNav(k),e.scales.xNav(C)]}))}},k.prototype._aaToNtCoord=function(t,e){return"-"===this.geneMap[t].strand?this.geneMap[t].end-3*e+1:this.geneMap[t].start+3*e-2},k.prototype._getZoomCoordinates=function(t,e){var n=[0,this.scales.xNav.domain()[1]],a=0;if(t.aa){var s=t.gene;a=1*((n=[e[s].start,e[s].end])[1]-n[0])}else{var o=this.scales.xNav.domain()[1];if(1===t.positions.length){var r=t.positions[0],i=.05*o;n=r-i<=0?[0,r+i]:r+i>=o?[r-i,o]:[r-i,r+i]}else{var c=Math.min.apply(null,t.positions),l=Math.max.apply(null,t.positions);n=[c-.05*(l-c),l+.05*(l-c)]}}return[Math.max(n[0]-a,0),Math.min(n[1]+a,this.scales.xNav.domain()[1])]},k.prototype._getSelectedNodes=function(t){if(this.aa===t.aa){var e=[];if(this.aa){var n,a=[],s=O(t.positions);try{for(s.s();!(n=s.n()).done;){var o=n.value;a.push([t.gene,o])}}catch(t){s.e(t)}finally{s.f()}var r,i=O(this.bars);try{for(i.s();!(r=i.n()).done;){var c,l=r.value,h=O(a);try{for(h.s();!(c=h.n()).done;){var p=c.value;l.prot===p[0]&&l.codon===p[1]&&e.push(l)}}catch(t){h.e(t)}finally{h.f()}}}catch(t){i.e(t)}finally{i.f()}}else{var u,d=O(this.bars);try{for(d.s();!(u=d.n()).done;){var f=u.value;-1!==t.positions.indexOf(f.x)&&e.push(f)}}catch(t){d.e(t)}finally{d.f()}}return e}console.error("entropy out of sync")},k.prototype._drawZoomGenes=function(t){var e=this;this.geneGraph.selectAll("*").remove();var n=this.scales.xNav.domain()[1]-this.scales.xNav.domain()[0]<1e6?"white":"black",a=!!t.filter((function(t){return"-"===t.strand})).length,s=function(t){return a?"-"===t?20:0:10},o=t.filter((function(t){return t.start<e.scales.xMain.domain()[1]&&t.start>e.scales.xMain.domain()[0]||t.end>e.scales.xMain.domain()[0]&&t.end<e.scales.xMain.domain()[1]||t.start<=e.scales.xMain.domain()[0]&&t.end>=e.scales.xMain.domain()[1]})),r=function(t){return t.start>e.scales.xMain.domain()[0]?e.scales.xMain(t.start):e.offsets.x1},i=function(t){return t.end<e.scales.xMain.domain()[1]?e.scales.xMain(t.end):e.offsets.x2},c=this.geneGraph.selectAll(".gene").data(o).enter().append("g");c.append("rect").attr("class","gene").attr("x",(function(t){return r(t)})).attr("y",(function(t){return s(t.strand)})).attr("width",(function(t){return i(t)-r(t)})).attr("height",20).style("fill",(function(t){return t.fill})).style("stroke",(function(){return n})),c.append("text").attr("x",(function(t){return e.scales.xMain(t.start)+(e.scales.xMain(t.end)-e.scales.xMain(t.start))/2})).attr("y",(function(t){return s(t.strand)+5})).attr("dy",".7em").attr("text-anchor","middle").style("fill",(function(){return"white"})).text((function(t){return i(t)-r(t)>15?t.prot:""}))},k.prototype._drawGenes=function(t){var e=this,n=this.scales.xNav.domain()[1]-this.scales.xNav.domain()[0]<1e6?"white":"black",a=this.navGraph.selectAll(".gene").data(t).enter().append("g");a.append("rect").attr("class","gene").attr("x",(function(t){return e.scales.xNav(t.start)})).attr("y",(function(t){return t.strand,5})).attr("width",(function(t){return e.scales.xNav(t.end)-e.scales.xNav(t.start)})).attr("height",20).style("fill",(function(t){return t.fill})).style("stroke",(function(){return n})),a.append("text").attr("x",(function(t){return e.scales.xNav(t.start)+(e.scales.xNav(t.end)-e.scales.xNav(t.start))/2})).attr("y",(function(t){return t.strand,10})).attr("dy",".7em").attr("text-anchor","middle").style("fill",(function(){return"white"})).text((function(t){return function(t){return t.end<e.scales.xNav.domain()[1]?e.scales.xNav(t.end):e.offsets.x2}(t)-function(t){return t.start>e.scales.xNav.domain()[0]?e.scales.xNav(t.start):e.offsets.x1}(t)>10?t.prot:""}))},k.prototype._clearSelectedBars=function(){var t,e=this,n=O(this.selectedNodes);try{for(n.s();!(t=n.n()).done;){var a=t.value,s=this.aa?"#prot".concat(a.prot).concat(a.codon):"#nt".concat(a.x),o=this.aa?function(t){return e.geneMap[t.prot].idx%2?u.h:u.a}:function(t){return t.prot?e.geneMap[t.prot].idx%2?u.h:u.a:u.f};Object(c.a)(s).style("fill",o)}}catch(t){n.e(t)}finally{n.f()}this.selectedNodes=[]},k.prototype._highlightSelectedBars=function(){var t,e=O(this.selectedNodes);try{for(e.s();!(t=e.n()).done;){var n=t.value;if(this.aa&&!n.prot)return;var a=this.aa?"#prot".concat(n.prot).concat(n.codon):"#nt".concat(n.x),s=this.aa||n.prot?this.geneMap[n.prot].fill:"red";Object(c.a)(a).style("fill",s)}}catch(t){e.e(t)}finally{e.f()}},k.prototype._drawBars=function(){var t=this;if(this.okToDrawBars){this.mainGraph.selectAll("*").remove();var e,n=this.scales.xMain.domain()[1]-this.scales.xMain.domain()[0];this.aa&&(n/=3),e=this.aa?n>600?2:function(e){return t.scales.xMain(t._aaToNtCoord(e.prot,e.codon)+2.6)-t.scales.xMain(t._aaToNtCoord(e.prot,e.codon))}:n>1e3?2:n>250?3:function(e){return t.scales.xMain(e.x+.3)-t.scales.xMain(e.x-.3)};var a=this.mainGraph.append("g").attr("clip-path","url(#clip)").selectAll(".bar"),s=this.aa?function(t){return"prot"+t.prot+t.codon}:function(t){return"nt"+t.x},o=this.aa?function(e){return t.scales.xMain(t._aaToNtCoord(e.prot,e.codon)-.3)}:function(e){return t.scales.xMain(e.x-.3)},r=this.aa?function(e){return t.geneMap[e.prot].idx%2?u.h:u.a}:function(e){return e.prot?t.geneMap[e.prot].idx%2?u.h:u.a:u.f};a.data(this.bars).enter().append("rect").attr("class","bar").attr("id",s).attr("x",o).attr("y",(function(e){return t.scales.y(e.y)})).attr("width",e).attr("height",(function(e){return t.offsets.heightMain-t.scales.y(e.y)})).style("fill",r).on("mouseover",(function(e){t.callbacks.onHover(e,d.c.pageX,d.c.pageY)})).on("mouseout",(function(e){t.callbacks.onLeave(e)})).on("click",(function(e){t.callbacks.onClick(e)})).style("cursor","pointer"),this._highlightSelectedBars()}},k.prototype._setScales=function(t,e){this.scales={},this.scales.xMax=t,this.scales.yMax=e,this.scales.yMin=0,this.scales.xMin=0,this.scales.xMain=Object(f.a)().domain([0,t]).range([this.offsets.x1,this.offsets.x2]),this.scales.xNav=Object(f.a)().domain([0,t]).range([this.offsets.x1,this.offsets.x2]),this.scales.y=Object(f.a)().domain([this.scales.yMin,1.2*e]).range([this.offsets.y2Main,this.offsets.y1Main])},k.prototype._drawAxes=function(){this.axes={},this.axes.y=Object(y.b)(this.scales.y).ticks(4),this.axes.xMain=Object(y.a)(this.scales.xMain).ticks(16),this.axes.xNav=Object(y.a)(this.scales.xNav).ticks(20,",f"),this.scales.xNav.domain()[1]-this.scales.xNav.domain()[0]>1e6&&this.axes.xNav.tickFormat(Object(v.a)(".1e")),this.svg.append("g").attr("class","y axis").attr("id","entropyYAxis").attr("transform","translate("+(this.offsets.x1+15)+","+this.offsets.y1Main+")").call(this.axes.y),this.svg.append("g").attr("class","xMain axis").attr("transform","translate("+this.offsets.x1+","+this.offsets.y2Main+")").call(this.axes.xMain),this.svg.append("g").attr("class","xNav axis").attr("transform","translate("+this.offsets.x1+","+this.offsets.y2Nav+")").call(this.axes.xNav)},k.prototype._updateYScaleAndAxis=function(t){this.scales.y=Object(f.a)().domain([this.scales.yMin,1.2*t]).range([this.offsets.y2Main,this.offsets.y1Main]),this.axes.y=Object(y.b)(this.scales.y).ticks(4),this.svg.select("#entropyYAxis").remove(),this.svg.append("g").attr("class","y axis").attr("id","entropyYAxis").attr("transform","translate("+(this.offsets.x1+15)+","+this.offsets.y1Main+")").call(this.axes.y)},k.prototype._calcOffsets=function(t,e){this.offsets={x1:15,x2:t-32,y1Main:0,y1Nav:e-65,y2Main:e-130,y2Nav:e-35,y1Gene:e-107,y2Gene:e-95},this.offsets.heightMain=this.offsets.y2Main-this.offsets.y1Main,this.offsets.heightNav=this.offsets.y2Nav-this.offsets.y1Nav,this.offsets.heightGene=this.offsets.y2Gene-this.offsets.y1Gene,this.offsets.width=this.offsets.x2-this.offsets.x1},k.prototype._addBrush=function(){var t=this;this.brushed=function(){var t=this,e=(d.c.selection||this.scales.xNav.range()).map(this.scales.xNav.invert,this.scales.xNav);this.zoomCoordinates=e.map(Math.round),d.c.selection?this._zoom(e[0],e[1]):this.navGraph.select(".brush").call(this.brush.move,(function(){return t.zoomCoordinates=t.scales.xNav.range().map(Math.round),t.scales.xNav.range()}))},this.brushFinished=function(){var t;this.brushed(),d.c.sourceEvent instanceof MouseEvent?d.c.selection&&"d3entropyParent"!==d.c.sourceEvent.target.id&&""!==d.c.sourceEvent.target.id?(d.c.sourceEvent.target.id.match(/^prot/)||d.c.sourceEvent.target.id.match(/^nt/))&&this.props.dispatch(Object(w.b)([void 0,void 0])):this.props.dispatch(Object(w.b)(this.zoomCoordinates)):(t=d.c)&&t.sourceEvent&&"zoom"===t.sourceEvent.type&&this.props.dispatch(Object(w.b)(this.zoomCoordinates))},this._zoom=function(t,e){var n=this,a=[t,e];this.scales.xMain.domain(a),this.axes.xMain=this.axes.xMain.scale(this.scales.xMain),this.svg.select(".xMain.axis").call(this.axes.xMain),this._drawBars(),this._drawZoomGenes(this.annotations),this.brushHandle&&this.brushHandle.attr("display",null).attr("transform",(function(t,e){return"translate("+n.scales.xNav(a[e])+","+(n.offsets.heightNav+25)+")"}))},this.brush=Object(x.a)().extent([[this.offsets.x1,0],[this.offsets.width+20,this.offsets.heightNav-1+25]]).on("brush",(function(){t.brushed()})).on("end",(function(){t.brushFinished()})),this.gBrush=this.navGraph.append("g").attr("class","brush").attr("stroke-width",0).call(this.brush).call(this.brush.move,(function(){return t.zoomCoordinates.map(t.scales.xNav)})),this.brushHandle=this.gBrush.selectAll(".handle--custom").data([{type:"w"},{type:"e"}]).enter().append("path").attr("class","handle--custom").attr("fill",u.a).attr("cursor","ew-resize").attr("d","M0,0 0,0 -5,11 5,11 0,0 Z").attr("transform",(function(e){return"e"===e.type?"translate("+(t.scales.xNav(t.zoomCoordinates[1])-1)+","+(t.offsets.heightNav+25)+")":"translate("+(t.scales.xNav(t.zoomCoordinates[0])+1)+","+(t.offsets.heightNav+25)+")"}))},k.prototype._addZoomLayers=function(){var t=this,e=[[this.offsets.x1,this.offsets.y1Main],[this.offsets.width,this.offsets.y2Main]];this.zoom=Object(m.a)().translateExtent(e).extent(e).on("zoom",(function(){return t.zoomed()}));var n=["option","shift"];g.a.bind(n,(function(){t.svg.append("rect").attr("class","overlay").attr("text","zoom").attr("transform","translate("+t.offsets.x1+","+t.offsets.y1Main+")").attr("width",t.offsets.width).attr("height",t.offsets.y2Nav+30-t.offsets.y1Main).call(t.zoom).on("wheel",(function(){d.c.preventDefault()}))}),"keydown"),g.a.bind(n,(function(){t.svg.selectAll(".overlay").remove(),t.svg.selectAll(".brush").remove(),t._addBrush()}),"keyup")},k.prototype._createZoomFn=function(){return function(){var t=this,e=d.c.transform,n=this.zoomCoordinates[1]-this.zoomCoordinates[0],a=(n-n/e.k)/2,s=[Math.max(this.zoomCoordinates[0]+a,this.scales.xNav(0)),Math.min(this.zoomCoordinates[1]-a,this.scales.xNav.domain()[1])];s[1]-s[0]<500&&e.k>1||(this.zoomCoordinates=s,e.rescaleX(this.scales.xMain),this.axes.xMain=this.axes.xMain.scale(this.scales.xMain),this.svg.select(".xMain.axis").call(this.axes.xMain),this._drawBars(),this._drawZoomGenes(this.annotations),this.navGraph.select(".brush").call(this.brush.move,(function(){return t.zoomCoordinates.map(t.scales.xNav)})))}},k.prototype._drawMainNavElements=function(){this.mainGraph=this.svg.append("g").attr("class","main").attr("transform","translate("+this.offsets.x1+","+this.offsets.y1Main+")"),this.navGraph=this.svg.append("g").attr("class","nav").attr("transform","translate("+this.offsets.x1+","+this.offsets.y1Nav+")"),this.geneGraph=this.svg.append("g").attr("class","Gene").attr("transform","translate("+this.offsets.x1+","+this.offsets.y1Gene+")")},k.prototype._addClipMask=function(){this.svg.append("g").append("clipPath").attr("class","clipPath").attr("transform","translate("+this.offsets.x1+","+this.offsets.y1Main+")").attr("id","clip").append("rect").attr("id","cliprect").attr("width",this.offsets.width).attr("height",this.offsets.heightMain)};var C,_,z=k,E=function(t){var e=t.hovered,n=t.width,a=t.height,o=t.mutType,r=t.showCounts,i=t.geneMap,c=t.t;if(!e)return null;var l={container:{position:"absolute",width:200,padding:"10px",borderRadius:10,zIndex:1e3,pointerEvents:"none",fontFamily:u.e.panel.fontFamily,fontSize:14,lineHeight:1,fontWeight:u.e.panel.fontWeight,color:u.e.panel.color,backgroundColor:u.e.panel.backgroundColor,wordWrap:"break-word",wordBreak:"break-word"}},h=document.getElementById("entropyYAxis").getBoundingClientRect(),p={x:e.x-h.left+15,y:e.y-h.top};p.x<.7*n?l.container.left=p.x:l.container.right=n-p.x,p.y<.5*a?l.container.top=p.y:l.container.bottom=a-p.y;var d=e.d.prot?"-"===i[e.d.prot].strand:null,f=e.d.prot?"aa"===o?d?i[e.d.prot].end-3*e.d.codon+3:i[e.d.prot].start+3*e.d.codon:d?i[e.d.prot].end-e.d.x:e.d.x-i[e.d.prot].start-1:null,y=e.d.prot&&"aa"!==o?Math.floor(f/3)+1:null;return s.a.createElement("div",{style:l.container},s.a.createElement("div",{className:"tooltip",style:u.e.tooltip},s.a.createElement("div",null,"aa"===o?c("Codon {{codon}} in protein {{protein}}",{codon:e.d.codon,protein:e.d.prot}):e.d.prot?"".concat(c("Nucleotide {{nuc}}",{nuc:e.d.x})," (").concat(c("Codon {{codon}} in protein {{protein}}",{codon:y,protein:e.d.prot}),")"):c("Nucleotide {{nuc}}",{nuc:e.d.x})),s.a.createElement("p",null),s.a.createElement("div",null,"aa"===o?c("Nuc positions {{a}} to {{b}}",{a:f-2,b:f}):""),s.a.createElement("p",null),s.a.createElement("div",null,null===d?"":c(d?"Negative strand":"Positive strand")),s.a.createElement("p",null),s.a.createElement("div",null,r?"".concat(c("Num mutations"),": ").concat(e.d.y):"".concat(c("entropy"),": ").concat(e.d.y)),s.a.createElement("div",{style:u.e.comment},c("Click to color tree & map"))))},S=n(830);n(838),n(920);function B(t){return(B="function"==typeof Symbol&&"symbol"==typeof Symbol.iterator?function(t){return typeof t}:function(t){return t&&"function"==typeof Symbol&&t.constructor===Symbol&&t!==Symbol.prototype?"symbol":typeof t})(t)}function T(t,e){var n=Object.keys(t);if(Object.getOwnPropertySymbols){var a=Object.getOwnPropertySymbols(t);e&&(a=a.filter((function(e){return Object.getOwnPropertyDescriptor(t,e).enumerable}))),n.push.apply(n,a)}return n}function A(t){for(var e=1;e<arguments.length;e++){var n=null!=arguments[e]?arguments[e]:{};e%2?T(Object(n),!0).forEach((function(e){D(t,e,n[e])})):Object.getOwnPropertyDescriptors?Object.defineProperties(t,Object.getOwnPropertyDescriptors(n)):T(Object(n)).forEach((function(e){Object.defineProperty(t,e,Object.getOwnPropertyDescriptor(n,e))}))}return t}function G(t,e){for(var n=0;n<e.length;n++){var a=e[n];a.enumerable=a.enumerable||!1,a.configurable=!0,"value"in a&&(a.writable=!0),Object.defineProperty(t,Z(a.key),a)}}function P(t,e){return(P=Object.setPrototypeOf?Object.setPrototypeOf.bind():function(t,e){return t.__proto__=e,t})(t,e)}function L(t){var e=function(){if("undefined"==typeof Reflect||!Reflect.construct)return!1;if(Reflect.construct.sham)return!1;if("function"==typeof Proxy)return!0;try{return Boolean.prototype.valueOf.call(Reflect.construct(Boolean,[],(function(){}))),!0}catch(t){return!1}}();return function(){var n,a=Y(t);if(e){var s=Y(this).constructor;n=Reflect.construct(a,arguments,s)}else n=a.apply(this,arguments);return R(this,n)}}function R(t,e){if(e&&("object"===B(e)||"function"==typeof e))return e;if(void 0!==e)throw new TypeError("Derived constructors may only return object or undefined");return function(t){if(void 0===t)throw new ReferenceError("this hasn't been initialised - super() hasn't been called");return t}(t)}function Y(t){return(Y=Object.setPrototypeOf?Object.getPrototypeOf.bind():function(t){return t.__proto__||Object.getPrototypeOf(t)})(t)}function D(t,e,n){return(e=Z(e))in t?Object.defineProperty(t,e,{value:n,enumerable:!0,configurable:!0,writable:!0}):t[e]=n,t}function Z(t){var e=function(t,e){if("object"!==B(t)||null===t)return t;var n=t[Symbol.toPrimitive];if(void 0!==n){var a=n.call(t,e||"default");if("object"!==B(a))return a;throw new TypeError("@@toPrimitive must return a primitive value.")}return("string"===e?String:Number)(t)}(t,"string");return"symbol"===B(e)?e:String(e)}var H=Object(i.b)((function(t){return{mutType:t.controls.mutType,bars:t.entropy.bars,annotations:t.entropy.annotations,geneMap:t.entropy.geneMap,geneLength:t.controls.geneLength,maxYVal:t.entropy.maxYVal,showCounts:t.entropy.showCounts,loaded:t.entropy.loaded,colorBy:t.controls.colorBy,zoomMin:t.controls.zoomMin,zoomMax:t.controls.zoomMax,defaultColorBy:t.controls.defaults.colorBy,shouldReRender:!1,panelLayout:t.controls.panelLayout,narrativeMode:t.narrative.display}}))((D(_=function(t){!function(t,e){if("function"!=typeof e&&null!==e)throw new TypeError("Super expression must either be null or a function");t.prototype=Object.create(e&&e.prototype,{constructor:{value:t,writable:!0,configurable:!0}}),Object.defineProperty(t,"prototype",{writable:!1}),e&&P(t,e)}(r,t);var e,n,a,o=L(r);function r(t){var e;return function(t,e){if(!(t instanceof e))throw new TypeError("Cannot call a class as a function")}(this,r),(e=o.call(this,t)).state={hovered:!1,chart:!1},e}return e=r,(n=[{key:"onHover",value:function(t,e,n){this.setState({hovered:{d:t,type:".tip",x:e,y:n}})}},{key:"onLeave",value:function(){this.setState({hovered:!1})}},{key:"onClick",value:function(t){if(!this.props.narrativeMode){var e=function(t,e){return Object(M.d)("aa"===t?{gene:e.prot,positions:[e.codon]}:{positions:[e.x]})}(this.props.mutType,t);Object(S.a)("color-by-genotype"),this.props.dispatch(Object(p.a)(e)),this.setState({hovered:!1})}}},{key:"changeMutTypeCallback",value:function(t){t!==this.props.mutType&&(this.props.dispatch(Object(p.a)(this.props.defaultColorBy)),this.props.dispatch(Object(w.a)(t)))}},{key:"aaNtSwitch",value:function(t){var e=this;return this.props.narrativeMode?null:s.a.createElement("div",{style:A(A({},u.k),t.aaNtSwitch)},s.a.createElement("button",{key:1,style:"aa"===this.props.mutType?u.m:u.l,onClick:function(){return e.changeMutTypeCallback("aa")}},s.a.createElement("span",{style:t.switchTitle}," ","AA"," ")),s.a.createElement("button",{key:2,style:"aa"!==this.props.mutType?u.m:u.l,onClick:function(){return e.changeMutTypeCallback("nuc")}},s.a.createElement("span",{style:t.switchTitle}," ","NT"," ")))}},{key:"entropyCountSwitch",value:function(t){var e=this,n=this.props.t;return this.props.narrativeMode?null:s.a.createElement("div",{style:A(A({},u.k),t.entropyCountSwitch)},s.a.createElement("button",{key:1,style:this.props.showCounts?u.l:u.m,onClick:function(){return e.props.dispatch(Object(w.c)(!1))}},s.a.createElement("span",{style:t.switchTitle}," ",n("entropy")," ")),s.a.createElement("button",{key:2,style:this.props.showCounts?u.m:u.l,onClick:function(){return e.props.dispatch(Object(w.c)(!0))}},s.a.createElement("span",{style:t.switchTitle}," ",n("events")," ")))}},{key:"setUp",value:function(t){var e=new z(this.d3entropy,t.annotations,t.geneMap,t.geneLength[N.K],{onHover:this.onHover.bind(this),onLeave:this.onLeave.bind(this),onClick:this.onClick.bind(this)});e.render(t),t.narrativeMode&&Object(c.a)(this.d3entropy).selectAll(".handle--custom").style("visibility","hidden"),this.setState({chart:e})}},{key:"componentDidMount",value:function(){this.props.loaded&&this.setUp(this.props)}},{key:"componentWillReceiveProps",value:function(t){if(t.loaded||this.setState({chart:!1}),this.state.chart){if(this.props.width!==t.width||this.props.height!==t.height)this.state.chart.render(t);else{var e={};if(this.props.zoomMax===t.zoomMax&&this.props.zoomMin===t.zoomMin||(e.zoomMax=t.zoomMax,e.zoomMin=t.zoomMin),this.props.bars!==t.bars&&(e.aa="aa"===t.mutType,e.newBars=t.bars,e.maxYVal=t.maxYVal),this.props.colorBy!==t.colorBy&&(Object(M.f)(this.props.colorBy)||Object(M.f)(t.colorBy)))if(Object(M.f)(t.colorBy)){var n=Object(M.a)(t.colorBy,t.geneLength);if(n.aa)e.gene=n.gene,e.start=t.geneMap[e.gene].start,e.end=t.geneMap[e.gene].end;else{var a,s,o=n.positions,r=this.state.chart.zoomCoordinates,i=this.state.chart.maxNt;if(o.length>1){var l=Math.min.apply(null,o),h=Math.max.apply(null,o);a=l-.05*(h-l),s=h+.05*(h-l)}else{var p=o[0],u=.05*i,d=p-u<=0?[0,p+u]:p+u>=i?[p-u,i]:[p-u,p+u];a=d[0],s=d[1]}a>r[0]-.4*i&&a<r[0]+.4*i&&s>r[1]-.4*i&&s<r[1]+.4*i&&o.every((function(t){return t>r[0]}))&&o.every((function(t){return t<r[1]}))||(e.gene=n.gene,e.start=a,e.end=s)}e.selected=Object(M.a)(t.colorBy,t.geneLength)}else e.clearSelected=!0;Object.keys(e).length&&this.state.chart.update(e)}this.props.narrativeMode!==t.narrativeMode&&(t.narrativeMode?Object(c.a)(this.d3entropy).selectAll(".handle--custom").style("visibility","hidden"):Object(c.a)(this.d3entropy).selectAll(".handle--custom").style("visibility","visible"))}else t.loaded&&this.setUp(t)}},{key:"render",value:function(){var t,e=this,n=this.props.t,a={switchContainer:{position:"absolute",marginTop:-5,paddingLeft:(t=this.props.width)-100},switchContainerWide:{position:"absolute",marginTop:-25,paddingLeft:t-185},switchTitle:{margin:5,position:"relative",top:-1},aaNtSwitch:{position:"absolute",right:5,top:0,zIndex:100},entropyCountSwitch:{position:"absolute",right:74,top:0,zIndex:100}};return s.a.createElement(h.a,{title:n("Diversity")},s.a.createElement(E,{hovered:this.state.hovered,width:this.props.width,height:this.props.height,mutType:this.props.mutType,showCounts:this.props.showCounts,geneMap:this.props.geneMap,t:n}),s.a.createElement("svg",{id:"d3entropyParent",style:{pointerEvents:"auto"},width:this.props.width,height:this.props.height},s.a.createElement("g",{ref:function(t){e.d3entropy=t},id:"d3entropy"})),this.aaNtSwitch(a),this.entropyCountSwitch(a))}},{key:"componentWillUnmount",value:function(){}}])&&G(e.prototype,n),a&&G(e,a),Object.defineProperty(e,"prototype",{writable:!1}),r}(s.a.Component),"propTypes",{dispatch:r.a.func.isRequired,loaded:r.a.bool.isRequired,colorBy:r.a.string.isRequired,defaultColorBy:r.a.string.isRequired,mutType:r.a.string.isRequired}),C=_))||C,F=Object(l.a)()(H);e.default=F},920:function(t,e,n){var a=n(921);"string"==typeof a&&(a=[[t.i,a,""]]);n(235)(a,{});a.locals&&(t.exports=a.locals)},921:function(t,e,n){(t.exports=n(338)(!1)).push([t.i,'\n/* Variables - keep in sync with globalStyles.js*/\n:root {\n  --darkGrey: #333;\n  --medGrey: #888;\n}\n\n#d3entropy .overlay {\n  fill: none;\n  pointer-events: all;\n}\n\n#d3entropy .brush {\n  stroke: none;\n}\n\n#d3entropy .niceText {\n  font-family: "Lato", "Helvetica Neue", "Helvetica", "sans-serif";\n  font-size: 14px;\n  font-style: italic;\n  color: var(--medGrey);\n}\n\n.brush .selection {\n  cursor: grab;\n  cursor: -moz-grab;\n  cursor: -webkit-grab;\n}\n',""])}}]);